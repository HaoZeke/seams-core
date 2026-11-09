"""i-PI socket client. Positions on the wire are bohr."""

import socket

import numpy as np

BOHR_A = 0.529177210903
EV_PER_HA = 27.211386245988
HEADER = 12


def message(text):
    raw = text.encode("ascii")
    if len(raw) > HEADER:
        raise ValueError(text)
    return raw + b" " * (HEADER - len(raw))


def read_exact(sock, n):
    buf = bytearray()
    while len(buf) < n:
        chunk = sock.recv(n - len(buf))
        if not chunk:
            raise ConnectionError("i-PI socket closed")
        buf += chunk
    return bytes(buf)


def read_array(sock, dtype, count):
    raw = read_exact(sock, count * np.dtype(dtype).itemsize)
    return np.frombuffer(raw, dtype=dtype).copy()


def send_array(sock, array):
    sock.sendall(np.ascontiguousarray(array).tobytes())


def lj_forces(pos_a, cell_a, epsilon=0.01, sigma=2.5, cutoff=6.0):
    n = int(pos_a.shape[0])
    force = np.zeros_like(pos_a)
    energy = 0.0
    inv = np.linalg.inv(cell_a)
    cut2 = cutoff * cutoff
    for i in range(n):
        for j in range(i + 1, n):
            delta = pos_a[j] - pos_a[i]
            frac = inv @ delta
            frac -= np.round(frac)
            mic = cell_a @ frac
            r2 = float(mic @ mic)
            if r2 <= 0.0 or r2 > cut2:
                continue
            sr6 = (sigma * sigma / r2) ** 3
            energy += 4.0 * epsilon * (sr6 * sr6 - sr6)
            mag = 24.0 * epsilon * (2.0 * sr6 * sr6 - sr6) / r2
            force[i] -= mag * mic
            force[j] += mag * mic
    return energy, force


class IpiClient:
    def __init__(self, host, port):
        self.host = host
        self.port = int(port)
        self.initialised = False
        self.have_data = False
        self.nat = 0
        self.pending = None

    def run(self, on_frame, max_frames=0):
        sock = socket.create_connection((self.host, self.port))
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        frames = 0
        try:
            while max_frames <= 0 or frames < max_frames:
                header = read_exact(sock, HEADER)
                if header == message("STATUS"):
                    self._status(sock)
                elif header == message("INIT"):
                    self._init(sock)
                elif header == message("POSDATA"):
                    self._positions(sock, on_frame, frames)
                    frames += 1
                elif header == message("GETFORCE"):
                    self._forces(sock)
                elif header == message("EXIT"):
                    break
                else:
                    raise RuntimeError("unexpected i-PI header " + repr(header))
        finally:
            sock.close()
        return frames

    def _status(self, sock):
        if not self.initialised:
            sock.sendall(message("NEEDINIT"))
        elif self.have_data:
            sock.sendall(message("HAVEDATA"))
        else:
            sock.sendall(message("READY"))

    def _init(self, sock):
        read_array(sock, np.int32, 1)
        n = int(read_array(sock, np.int32, 1)[0])
        if n > 0:
            read_exact(sock, n)
        self.initialised = True

    def _positions(self, sock, on_frame, frames):
        cell = read_array(sock, np.float64, 9).reshape(3, 3)
        read_array(sock, np.float64, 9)
        self.nat = int(read_array(sock, np.int32, 1)[0])
        pos = read_array(sock, np.float64, self.nat * 3).reshape(self.nat, 3)
        cell_a = cell * BOHR_A
        pos_a = pos * BOHR_A
        on_frame(frames, pos_a, cell_a)
        energy_ev, force_ev_a = lj_forces(pos_a, cell_a)
        scale = BOHR_A / EV_PER_HA
        self.pending = (energy_ev / EV_PER_HA, force_ev_a * scale)
        self.have_data = True

    def _forces(self, sock):
        if self.pending is None:
            raise RuntimeError("GETFORCE before positions")
        energy, force = self.pending
        sock.sendall(message("FORCEREADY"))
        send_array(sock, np.float64(energy))
        send_array(sock, np.int32(self.nat))
        send_array(sock, np.ascontiguousarray(force, dtype=np.float64))
        send_array(sock, np.zeros((3, 3), dtype=np.float64))
        send_array(sock, np.int32(0))
        self.have_data = False
        self.pending = None
