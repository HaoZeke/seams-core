library(dseams)

root <- Sys.getenv("DSEAMS_ROOT", unset = "")
tip <- file.path(root, "input", "traj", "exampleTraj.lammpstrj")
cubic <- file.path(root, "input", "traj", "mW_cubic.lammpstrj")

skip_if(root == "" || !file.exists(tip), "DSEAMS_ROOT is not set")

test_that("CHILL+ on the mixed TIP4P frame matches the CLI", {
  frame <- read_frame(tip, type = 2L)
  expect_equal(nrow(positions(frame)), 250L)
  ice <- chill_plus(frame)
  expect_s3_class(ice, "data.frame")
  tab <- table(ice$name)
  expect_false("clathrate" %in% names(tab))
  expect_equal(as.integer(tab[["interClathrate"]]), 9L)
  expect_equal(as.integer(tab[["water"]]), 241L)
  expect_identical(ice$label, chill_plus(frame, graph = "cutoff")$label)
  mutual <- table(chill_plus(frame, graph = "mutual")$name)
  expect_equal(as.integer(mutual[["clathrate"]]), 1L)
  expect_equal(as.integer(mutual[["interClathrate"]]), 7L)
  expect_equal(as.integer(mutual[["water"]]), 242L)
  cages_df <- cages(frame)
  expect_equal(attr(cages_df, "n_water"), 250L)
  expect_equal(attr(cages_df, "n_ih"), 0L)
  rings <- primitive_rings(frame, max_depth = 8L)
  expect_equal(rings$n[rings$size == 6L], 224L)
  sig <- signature_cages(frame, "512")
  expect_true(all(sig >= 0L))
  expect_length(symmetry(frame), 250L)
  q3 <- steinhardt(frame, order_l = 3L)
  # The default graph is the cutoff shell (seams_frame_steinhardt).
  # 0.389362 is the mutual four-nearest mean and must not be this value.
  got <- mean(q3$ql[is.finite(q3$ql)])
  expect_equal(got, 0.374897, tolerance = 1e-4)
  expect_false(abs(got - 0.389362) < 1e-4)
})

test_that("seeded cages under SEAMS_OFFLOAD match the host and the CLI", {
  nuc <- file.path(root, "input", "traj", "nucleation_f1.lammpstrj")
  skip_if(!file.exists(nuc))
  frame <- read_frame(nuc, type = 1L)
  Sys.unsetenv("SEAMS_OFFLOAD")
  host <- cages(frame, graph = "seeded")
  Sys.setenv(SEAMS_OFFLOAD = "1")
  off <- cages(frame, graph = "seeded")
  Sys.unsetenv("SEAMS_OFFLOAD")
  expect_identical(off$cage, host$cage)
  expect_equal(attr(off, "n_ih"), attr(host, "n_ih"))
  expect_equal(attr(off, "n_ic"), attr(host, "n_ic"))
  expect_equal(attr(off, "n_both"), attr(host, "n_both"))
  expect_equal(attr(off, "n_water"), attr(host, "n_water"))
  n <- attr(host, "n_ih") + attr(host, "n_ic") + attr(host, "n_both") +
    attr(host, "n_water")
  expect_equal(n, nrow(host))
  expect_gt(attr(host, "n_both"), 0L)
  bin <- Sys.getenv("SEAMS_BIN", unset = "")
  skip_if(bin == "" || !file.exists(bin))
  # R prepends its own lib directories to LD_LIBRARY_PATH. The engine is
  # linked against the toolchain libc, so those prefixes are removed for
  # this one call.
  saved_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
  parts <- strsplit(saved_ld, ":", fixed = TRUE)[[1]]
  parts <- parts[!grepl("^/usr/lib", parts)]
  Sys.setenv(LD_LIBRARY_PATH = paste(parts, collapse = ":"))
  line <- system2(bin, c("cages", nuc, "--graph", "seeded", "--cutoff", "3.5"),
                  stdout = TRUE)
  Sys.setenv(LD_LIBRARY_PATH = saved_ld)
  parts <- strsplit(line[[length(line)]], " +")[[1]]
  # nop N graph seeded hexagonal H cubic C both B water W
  expect_equal(as.integer(parts[[6]]), attr(host, "n_ih"))
  expect_equal(as.integer(parts[[8]]), attr(host, "n_ic"))
  expect_equal(as.integer(parts[[10]]), attr(host, "n_both"))
  expect_equal(as.integer(parts[[12]]), attr(host, "n_water"))
})

test_that("cubic ice is 4096 cubic and 4096 double-diamond", {
  skip_if(!file.exists(cubic))
  frame <- read_frame(cubic, type = 1L)
  ice <- chill_plus(frame)
  expect_equal(sum(ice$name == "cubic"), 4096L)
  found <- cages(frame)
  expect_equal(attr(found, "n_ic"), 4096L)
  expect_gte(nframes(cubic), 1L)
  expect_output(print(frame), "dseams frame")
  order <- q4_q6(frame)
  expect_equal(nrow(order), 4096L)
  tau <- jump_tau90(frame, frame, dt = 1)
  expect_true(is.na(tau) || is.finite(tau))
  expect_equal(nrow(steinhardt(frame, order_l = 6L, graph = "cutoff")), 4096L)
})

chordless_reference <- function(path) {
  raw <- scan(path, what = integer(), quiet = TRUE)
  atoms <- raw[[1]]
  n <- raw[[2]]
  size <- raw[[3]]
  body <- raw[-(1:3)]
  rows <- matrix(body[seq_len(n * size)], nrow = n, ncol = size, byrow = TRUE)
  list(atoms = atoms, size = size, rings = rows,
       per_atom = body[n * size + seq_len(atoms)])
}

test_that("rings_chordless matches the C reference", {
  ref_path <- Sys.getenv("SEAMS_CHORDLESS_REF", unset = "")
  traj <- Sys.getenv("SEAMS_CHORDLESS_TRAJ", unset = "")
  skip_if(ref_path == "" || !file.exists(ref_path), "SEAMS_CHORDLESS_REF is not set")
  ref <- chordless_reference(ref_path)
  frame <- read_frame(traj, frame = 1L, type = 2L)
  got <- rings_chordless(frame, size = ref$size, cutoff = 3.5, type = 0L)
  expect_identical(dim(got), dim(ref$rings))
  expect_identical(unname(got[, , drop = FALSE]) - 1L, ref$rings)
  expect_identical(attr(got, "per_atom"), ref$per_atom)
  expect_identical(positions(frame)[got[1, ], , drop = FALSE],
                   positions(frame)[ref$rings[1, ] + 1L, , drop = FALSE])
  sys <- seams_system_new()
  seams_system_read(sys, traj, frame = 1L, type = 2L)
  raw <- seams_rings_chordless(sys, cutoff = 3.5, size = ref$size)
  expect_identical(raw$rings, ref$rings)
  expect_identical(raw$per_atom, ref$per_atom)
})

test_that("cubic ice has twelve chordless six-cycles through each atom", {
  skip_if(!file.exists(cubic))
  frame <- read_frame(cubic, type = 1L)
  got <- rings_chordless(frame, size = 6L)
  expect_equal(nrow(got), 8192L)
  expect_true(all(attr(got, "per_atom") == 12L))
  expect_true(all(got[, 1] < apply(got[, -1], 1, min)))
  expect_equal(nrow(rings_chordless(frame, size = 4L)), 0L)
  expect_error(rings_chordless(frame, size = 13L))
})
