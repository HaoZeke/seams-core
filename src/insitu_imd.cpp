#include <seams_insitu.h>

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/time.h>
#include <vector>

namespace {

constexpr int kHeader = 8;
constexpr int kHandshake = 4;
constexpr int kEnergies = 1;
constexpr int kCoords = 2;
constexpr int kGo = 3;
constexpr int kVersion = 2;
constexpr int kEnergyBytes = 40;

bool read_full(int fd, void *data, int bytes) {
  auto *out = static_cast<unsigned char *>(data);
  int got = 0;
  while (got < bytes) {
    const ssize_t n = ::recv(fd, out + got, static_cast<std::size_t>(bytes - got), 0);
    if (n <= 0) {
      return false;
    }
    got += static_cast<int>(n);
  }
  return true;
}

bool write_full(int fd, const void *data, int bytes) {
  const auto *in = static_cast<const unsigned char *>(data);
  int sent = 0;
  while (sent < bytes) {
    const ssize_t n = ::send(fd, in + sent, static_cast<std::size_t>(bytes - sent), 0);
    if (n <= 0) return false;
    sent += static_cast<int>(n);
  }
  return true;
}

std::uint32_t from_be(std::uint32_t value) { return ntohl(value); }

} // namespace

int seams_insitu_imd_session(const char *host, int port, const char *dump, const char *record,
                             const double *box9, int max_frames) {
  if (host == nullptr || dump == nullptr || record == nullptr || box9 == nullptr || port <= 0) {
    return SEAMS_INSITU_ERR;
  }
  const int fd = ::socket(AF_INET, SOCK_STREAM, 0);
  timeval tv{};
  tv.tv_sec = 30;
  ::setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
  ::setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
  sockaddr_in addr{};
  addr.sin_family = AF_INET;
  addr.sin_port = htons(static_cast<uint16_t>(port));
  if (::inet_pton(AF_INET, host, &addr.sin_addr) != 1) {
    ::close(fd);
    return SEAMS_INSITU_ERR;
  }
  if (::connect(fd, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) != 0) {
    ::close(fd);
    return SEAMS_INSITU_ERR_IO;
  }
  unsigned char header[kHeader];
  if (!read_full(fd, header, kHeader)) {
    ::close(fd);
    return SEAMS_INSITU_ERR_IO;
  }
  std::uint32_t type_be = 0;
  std::uint32_t version = 0;
  std::memcpy(&type_be, header, 4);
  std::memcpy(&version, header + 4, 4);
  if (static_cast<int>(from_be(type_be)) != kHandshake || static_cast<int>(version) != kVersion) {
    ::close(fd);
    return SEAMS_INSITU_ERR;
  }
  std::uint32_t go[2] = {htonl(static_cast<std::uint32_t>(kGo)), htonl(0)};
  write_full(fd, go, kHeader);
  seams_insitu_engine *eng = seams_insitu_create(nullptr);
  if (seams_insitu_open_record(eng, record) != SEAMS_INSITU_OK) {
    seams_insitu_destroy(eng);
    ::close(fd);
    return SEAMS_INSITU_ERR_IO;
  }
  int frames = 0;
  int status = SEAMS_INSITU_OK;
  std::int32_t md_step = 0;
  while (max_frames <= 0 || frames < max_frames) {
    if (!read_full(fd, header, kHeader)) {
      break;
    }
    std::uint32_t raw_type = 0;
    std::uint32_t raw_len = 0;
    std::memcpy(&raw_type, header, 4);
    std::memcpy(&raw_len, header + 4, 4);
    const int kind = static_cast<int>(from_be(raw_type));
    const int length = static_cast<int>(from_be(raw_len));
    if (kind == kEnergies) {
      unsigned char body[kEnergyBytes];
      if (!read_full(fd, body, kEnergyBytes)) {
        status = SEAMS_INSITU_ERR_IO;
        break;
      }
      std::memcpy(&md_step, body, sizeof(md_step));
      continue;
    }
    if (kind != kCoords) {
      status = SEAMS_INSITU_ERR;
      break;
    }
    if (length < 0) {
      status = SEAMS_INSITU_ERR;
      break;
    }
    std::vector<float> packed(static_cast<std::size_t>(length) * 3);
    if (length > 0 &&
        !read_full(fd, packed.data(), static_cast<int>(packed.size() * sizeof(float)))) {
      status = SEAMS_INSITU_ERR_IO;
      break;
    }
    std::vector<double> xyz(packed.size());
    for (std::size_t i = 0; i < packed.size(); ++i) {
      xyz[i] = static_cast<double>(packed[i]);
    }
    std::int32_t step = md_step;
    seams_insitu_frame frame{};
    frame.step = step;
    frame.n = length;
    frame.xyz = xyz.empty() ? nullptr : xyz.data();
    frame.xlo = box9[0];
    frame.ylo = box9[1];
    frame.zlo = box9[2];
    frame.lx = box9[3];
    frame.ly = box9[4];
    frame.lz = box9[5];
    frame.xy = box9[6];
    frame.xz = box9[7];
    frame.yz = box9[8];
    frame.gathered = 1;
    const int rc = seams_insitu_observe(eng, &frame, nullptr);
    (void)dump;
    if (rc != SEAMS_INSITU_OK) {
      status = rc;
      break;
    }
    frames += 1;
  }
  seams_insitu_destroy(eng);
  ::close(fd);
  if (frames < 1 && status == SEAMS_INSITU_OK) {
    return SEAMS_INSITU_ERR_IO;
  }
  return status;
}

int seams_insitu_imd_main(int argc, char **argv) {
  if (argc != 15) {
    std::fprintf(stderr,
                 "usage: seams-gromacs-imd HOST PORT DUMP RECORD XLO YLO ZLO LX LY LZ XY XZ YZ MAX\n");
    return 2;
  }
  char *end = nullptr;
  const int port = static_cast<int>(std::strtol(argv[2], &end, 10));
  if (end == argv[2] || *end != '\0') {
    return 2;
  }
  double box[9];
  for (int i = 0; i < 9; ++i) {
    box[i] = std::strtod(argv[5 + i], &end);
    if (end == argv[5 + i] || *end != '\0') {
      return 2;
    }
  }
  const int max_frames = static_cast<int>(std::strtol(argv[14], &end, 10));
  if (end == argv[14] || *end != '\0') {
    return 2;
  }
  const int rc =
      seams_insitu_imd_session(argv[1], port, argv[3], argv[4], box, max_frames);
  return rc == SEAMS_INSITU_OK ? 0 : 1;
}
