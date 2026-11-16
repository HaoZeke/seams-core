# Read the GenIce sI cell: neighbours via steinhardt, then q6, CHILL+, rings,
# and chordless five-cycles.
using Dseams
using Statistics

path = ARGS[1]
frame = read_frame(path)
@assert nop(frame) == 46
order = steinhardt(frame, 6; cutoff = 3.5)
@assert length(order.ql) == 46
@assert mean(order.ql) > 0.5
ice = chill_plus(frame)
@assert ice.counts["clathrate"] == 46
rings = primitive_rings(frame; graph = "cutoff")
@assert rings[5] == 48
@assert rings[6] == 54
pentagons = rings_chordless(frame; size = 5)
@assert size(pentagons.cycles) == (5, 48)
@assert sum(pentagons.per_atom) == 5 * 48
println("julia example ok")
