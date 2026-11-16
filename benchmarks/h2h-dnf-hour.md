# Hour budget for the former 600 s DNF rows

Each call was allowed 3600 s. Peak RSS is the process high-water mark in kilobytes. Cores are 8. No host name.

| Kind | Molecules | Threads | Outcome | n | Other median ms | Peak RSS kb | What stopped it |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cycless | 884736 | 1 | measured | 8 | 2785794.3 | 6998288 | finished |
| cycless | 884736 | 8 | measured | 8 | 2774867.2 | 7349756 | finished |
| ira | 110592 | 1 | measured | 8 | 825641.7 | 174312 | finished |
| ira | 110592 | 8 | measured | 8 | 826970.3 | 158704 | finished |
| ira | 884736 | 1 | stopped | 0 |  | 320960 | exceeded 3600 s |
| ira | 884736 | 8 | stopped | 0 |  | 323892 | exceeded 3600 s. classify_serial has no thread argument, so the process was pinned to core 3 while the one-thread sessions that return inside the budget kept cores 0, 1 and 2 |
| sqq | 110592 | 1 | stopped | 0 |  | 422240 | exceeded 3600 s inside the cpp engine, which holds the interpreter so the in-process alarm could not fire; killed from outside. SQQ 0.5.8 engine cpp, ring.definition chordless. The python shortest_path engine does not return on these cells |
| sqq | 110592 | 8 | stopped | 0 |  | 423508 | exceeded 3600 s. The cpp extension does not start workers and holds the interpreter, so the outside killer stopped the process. It was pinned to core 4 while the one-thread sessions that return inside the budget kept cores 0, 1 and 2. SQQ 0.5.8 engine cpp, ring.definition chordless. The python shortest_path engine does not return on these cells |
| sqq | 4096 | 1 | measured | 8 | 1035697.9 | 229516 | finished |
| sqq | 4096 | 8 | measured | 8 | 1036923.2 | 235696 | finished |
| sqq | 884736 | 1 | stopped | 0 |  | 1915216 | exceeded 3600 s inside the cpp engine, which holds the interpreter so the in-process alarm could not fire; killed from outside. SQQ 0.5.8 engine cpp, ring.definition chordless. The python shortest_path engine does not return on these cells |
| sqq | 884736 | 8 | stopped | 0 |  | 1916176 | exceeded 3600 s. The cpp extension does not start workers and holds the interpreter, so the outside killer stopped the process. It was pinned to core 5 while the one-thread sessions that return inside the budget kept cores 0, 1 and 2. SQQ 0.5.8 engine cpp, ring.definition chordless. The python shortest_path engine does not return on these cells |

