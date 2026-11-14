# IRA and cycless scaling, labelled as fits

One sitting, core 0, one warm-up and three calls. IRA is classify_serial. cycless is cycles_iter on the 3.5 A graph, the graph built outside the timer, maximum size 6. Neither call takes a thread count, so the eight-thread DNF rows use this same serial fit.

The model is t = a N^b, ordinary least squares on log t against log N. a is exp(intercept). The standard error is the usual slope error. An extrapolated time is a fit, not a measurement. The 600 s budget still stands: a size that did not return is DNF.

| Kind | Molecules | Role | Median ms | Exponent | SE | R2 |
| --- | --- | --- | --- | --- | --- | --- |
| ira | 64 | measured | 11.6231 |  |  |  |
| cycless | 64 | measured | 157.7023 |  |  |  |
| ira | 216 | measured | 32.1538 |  |  |  |
| cycless | 216 | measured | 631.4641 |  |  |  |
| ira | 512 | measured | 81.3093 |  |  |  |
| cycless | 512 | measured | 1601.2396 |  |  |  |
| ira | 1000 | measured | 181.0865 |  |  |  |
| cycless | 1000 | measured | 3107.0783 |  |  |  |
| ira | 1728 | measured | 355.6018 |  |  |  |
| cycless | 1728 | measured | 5268.8832 |  |  |  |
| ira | 4096 | measured | 1236.5377 |  |  |  |
| cycless | 4096 | measured | 12380.6008 |  |  |  |
| ira | 8000 | measured | 3639.1208 |  |  |  |
| cycless | 8000 | measured | 24058.3087 |  |  |  |
| ira | 13824 | measured | 9727.2569 |  |  |  |
| cycless | 13824 | measured | 41399.2187 |  |  |  |
| ira | 32768 | measured | 51403.1226 |  |  |  |
| cycless | 32768 | measured | 98174.6126 |  |  |  |
| ira | 110592 | fit | 501919.9 | 1.8826 | 0.0354 | 0.999646 |
| ira | 884736 | fit | 25162336.6 | 1.8826 | 0.0354 | 0.999646 |
| ira | 110592 | quadratic | 585513.7 | 2 |  |  |
| ira | 884736 | quadratic | 37472876.4 | 2 |  |  |
| cycless | 110592 | fit | 354484.3 | 1.0193 | 0.0111 | 0.999166 |
| cycless | 884736 | fit | 2952280.4 | 1.0193 | 0.0111 | 0.999166 |

ira: b = 1.8826, se 0.0354, R2 0.999646, 3 sizes.

cycless: b = 1.0193, se 0.0111, R2 0.999166, 9 sizes.

cycless is close to linear on the whole ladder. The nine-point fit at 110592 molecules is 354484 ms. The summary sitting measured 341558 ms there, one thread, eight samples. The fit is not a replacement for that row. At 884736 molecules the same fit is 2952 s, past the 600 s budget, which is why that row is DNF.

IRA is not one exponent. cluster() subtracts every molecule and argsorts that vector, once per atom, so the leading cost is N^2. The local exponent rises from 0.84 between 64 and 216 molecules to 1.93 between 13824 and 32768. A nine-point line gives b = 1.346 and would put 110592 molecules near 142 s. That line is not quoted. The quoted fit is ordinary least squares on the three largest finished sizes, 8000, 13824, and 32768: b = 1.883, standard error 0.035. It puts classify_serial at 502 s for 110592 molecules and 25162 s for 884736. An N^2 projection of the 32768-molecule measurement, 51.403 s, is 586 s and 37473 s. The 600 s probe times the whole process and did not return. 502 s of classify_serial, still on a rising local exponent, sits next to that budget. Both numbers are fits, not measurements.
