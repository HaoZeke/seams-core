;; Copyright (c) d-SEAMS developers
;; Chordless six-cycles on the oxygen network of the mixed TIP4P frame.
(local dseams (require :dseams))
(local root (or (. arg 1) "."))
(local cloud (dseams.read (.. root "/input/traj/exampleTraj.lammpstrj") {:type 2}))
(local six (dseams.rings-chordless cloud {:size 6 :cutoff 3.5}))
(var through 0)
(each [_ n (ipairs six.per_atom)]
  (set through (+ through n)))
(assert (= through (* 6 six.n)) "each cycle passes six atoms")
(print (.. "tip chordless " six.n))
(cloud:close)
