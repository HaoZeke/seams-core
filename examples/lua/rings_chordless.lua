-- Copyright (c) d-SEAMS developers
-- Chordless five- and six-cycles of the GenIce sI cell.
local dseams = require("dseams")

local root = arg[1] or "."
local cloud = dseams.read(root .. "/input/traj/genice_sI.lammpstrj")
local five = dseams.rings_chordless(cloud, {size = 5, cutoff = 3.5})
local six = dseams.rings_chordless(cloud, {size = 6})
assert(five.n == 48 and six.n == 54, "sI has 48 pentagons and 54 hexagons")
local through = 0
for _, n in ipairs(five.per_atom) do
  through = through + n
end
assert(through == 5 * five.n)
local first = five.rings[1]
for k = 2, five.size do
  assert(first[k] > first[1], "each ring starts at its lowest atom")
end
print(string.format("chordless 5 %d 6 %d first %s", five.n, six.n, table.concat(first, " ")))
cloud:close()
