-- Copyright (c) d-SEAMS developers
-- LuaJIT FFI module generated from include/seams.h.
local ffi = require("ffi")
local C = dofile(arg[1])

local a = 2.75 * 4.0 / math.sqrt(3.0)
local fcc = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}
local basis = {{0, 0, 0}, {0.25, 0.25, 0.25}}
local xyz = ffi.new("double[24]")
local n = 0
for _, site in ipairs(fcc) do
  for _, shift in ipairs(basis) do
    xyz[3 * n] = (site[1] + shift[1]) * a
    xyz[3 * n + 1] = (site[2] + shift[2]) * a
    xyz[3 * n + 2] = (site[3] + shift[3]) * a
    n = n + 1
  end
end

local sys = C.seams_system_create()
assert(sys ~= nil)
local status = C.seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, nil, nil)
assert(status == 0, ffi.string(C.seams_last_error()))
local labels = ffi.new("int[8]")
local counts = ffi.new("int[9]")
status = C.seams_chill_plus(sys, 3.5, 0, 0, 0, labels, counts)
assert(status == 0, ffi.string(C.seams_last_error()))
assert(counts[0] == 8)
assert(C.seams_abi() == 1)
local cn = ffi.new("double[1]")
status = C.seams_coordination(sys, 3.5, 8, 0, 0, cn)
assert(status == 0, ffi.string(C.seams_last_error()))
local n_bonds = ffi.new("int[1]")
status = C.seams_hbonds(sys, 3.5, 1, 2, 2.42, 30, n_bonds)
assert(status == 0, ffi.string(C.seams_last_error()))
C.seams_system_destroy(sys)

local ref_path = os.getenv("SEAMS_CHORDLESS_REF")
local traj = os.getenv("SEAMS_CHORDLESS_TRAJ")
if ref_path ~= nil and ref_path ~= "" and traj ~= nil and traj ~= "" then
  local handle = assert(io.open(ref_path, "r"))
  local values = {}
  for token in handle:read("*a"):gmatch("%-?%d+") do
    values[#values + 1] = tonumber(token)
  end
  handle:close()
  local atoms, count, width = values[1], values[2], values[3]
  local water = C.seams_system_create()
  status = C.seams_system_read(water, traj, 1, 2)
  assert(status == 0, ffi.string(C.seams_last_error()))
  local nb = ffi.new("seams_neighbours *[1]")
  status = C.seams_neighbours_build(water, 3.5, 0, 0, nb)
  assert(status == 0, ffi.string(C.seams_last_error()))
  local got = ffi.new("size_t[1]")
  local per = ffi.new("int[?]", atoms)
  local rows = ffi.new("int[?]", math.max(count * width, 1))
  status = C.seams_rings_chordless(nb[0], width, rows, count * width, got, per, atoms)
  assert(status == 0, ffi.string(C.seams_last_error()))
  assert(tonumber(got[0]) == count, "ffi chordless count differs from the C reference")
  for i = 0, count * width - 1 do
    assert(rows[i] == values[4 + i], "ffi chordless row differs from the C reference")
  end
  for i = 0, atoms - 1 do
    assert(per[i] == values[4 + count * width + i], "ffi chordless per-atom differs")
  end
  C.seams_neighbours_destroy(nb[0])
  C.seams_system_destroy(water)
end
print("ffi ok")
