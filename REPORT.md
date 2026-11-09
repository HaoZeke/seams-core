# In-situ overhead

LAMMPS samples are its reported loop time divided by the step count.
Every other code times the same driver call and divides that wall clock
by how many steps ran.
The plain run and the in-situ run use the same input on one process and
one thread.
Six samples are kept.
The quartiles are Hyndman-Fan type 7.
The interquartile range (IQR) is the third quartile minus the first.
Overhead is the live sample minus the plain sample at the same index, in
seconds per step.

CPU `Intel(R) Xeon(R) Processor`.
Measured commit `ac1a97612f170ef03341e09084e38d0e049abf4a`.
Measured tree `61af403add2852df0ae2274ee8c0cea1dcb4ef41`.
libyodaLib SHA-256
`675190fd08612b2fc4da59085d3a66abcad3459f4b09476e31b74cfd92360f7e`.
LAMMPS SHA-256
`2ec47774f5ecd620929f4076ca033bbffec41a66f07331099ab3d0009f4cc5ad`.
GROMACS SHA-256
`dac429749a5a1dd5eb7c7c9249b59d9fd88cd454f46a959f7b5be3f53690020a`.
Plugin SHA-256
`1f6897cd7b105289877422ee24352d62ac70ea9638c881ed5179c07319ff7f93`.

| Code | Steps | Frames | Q1 overhead | Median | Q3 | IQR | Median ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lammps | 40 | 40 | 0.0014908149999999999 | 0.0015129174999999999 | 0.0015563537499999999 | 6.5538750000000007e-05 | 4.2962375975205358 |
| gromacs | 40 | 40 | 0.0013000000000000002 | 0.0013000000000000002 | 0.00131875 | 1.8749999999999887e-05 | 9.6666666666666661 |

Overhead is seconds per step.
The ratio is the live step time divided by the plain step time.
The GROMACS row uses the samples recorded for that code.
The size table records the native d-SEAMS module.
OpenMM writes a disk coordinate dump (DCD) in float32.
On this water system that file matches the in-situ record.

# LAMMPS size and cadence

The system is a diamond lattice of one site per water.
Plain and live samples are the LAMMPS loop time divided by the number of
steps.
The postprocess column is the wall time of the seams rings, chill-plus
and cages commands.
That column uses one native dump frame, repeated six times.
The in-situ run does not write a dump.
The byte count is the native dump of one frame.

| Waters | Cadence | Steps | Median absolute | Median ratio | Per analysed frame | Median postprocess | Dump bytes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1000 | 1 | 40 | 0.0015129174999999999 | 4.2962375975205358 | 0.0015129174999999999 | 0.024942176998592913 | 114451 |
| 1000 | 10 | 40 | 0.00044143124999999992 | 1.9712460363126132 | 0.004414312499999999 | 0.023702206002781168 | 114451 |
| 1000 | 100 | 100 | 0.00015351149999999995 | 1.3367028944178352 | 0.015351149999999996 | 0.02680031900308677 | 114451 |
| 64000 | 1 | 8 | 0.14838237500000001 | 9.6501131000035052 | 0.14838237500000001 | 0.50551247750263428 | 7546706 |
| 64000 | 10 | 10 | 0.070805950000000006 | 5.1706904309966291 | 0.70805950000000006 | 0.52629454850466573 | 7546706 |
| 64000 | 100 | 100 | 0.0082464500000000024 | 1.348126955754164 | 0.82464500000000029 | 0.51626300800126046 | 7546706 |

## lammps 1000 waters every 1

Warmup plain seconds per step:

0.00045032249999999997

Warmup live seconds per step:

0.0019604799999999997

Plain seconds per step:

0.00046736750000000006 0.00045604749999999996 0.00046040999999999996
0.00045836 0.00048961999999999996 0.00045961999999999999

Live seconds per step:

0.0020771775000000001 0.0020220849999999999 0.0019486525000000001
0.0019856624999999998 0.0019509925000000001 0.0019581525000000001

Postprocess wall seconds:

0.030568721005693078 0.025573327002348378 0.024311026994837448
0.02381401199818356 0.023269331999472342 0.025722330996359233

## lammps 1000 waters every 10

Warmup plain seconds per step:

0.00045732999999999999

Warmup live seconds per step:

0.00091365999999999995

Plain seconds per step:

0.00045482999999999998 0.00046161249999999998 0.000455795
0.00045321249999999999 0.00045577500000000001 0.0004581125

Live seconds per step:

0.0009478275 0.00089764500000000008 0.00089717249999999996
0.00089469749999999992 0.00090057750000000002 0.00089934750000000006

Postprocess wall seconds:

0.029020085996307898 0.022193682001670823 0.024095237000437919
0.023966407999978401 0.023438004005583934 0.023227479003253393

## lammps 1000 waters every 100

Warmup plain seconds per step:

0.00045477099999999997

Warmup live seconds per step:

0.00060882500000000001

Plain seconds per step:

0.000456687 0.00045519099999999997 0.000476174 0.00045485799999999999
0.00056981900000000005 0.00045453900000000001

Live seconds per step:

0.00060771899999999995 0.00061118199999999994 0.00075845200000000006
0.00061254799999999998 0.00061639699999999997 0.00060385899999999999

Postprocess wall seconds:

0.02999862300202949 0.026807374000782147 0.023080463994119782
0.024756757004070096 0.029882750000979286 0.026793264005391393

## lammps 64000 waters every 1

Warmup plain seconds per step:

0.017096375

Warmup live seconds per step:

0.18129124999999999

Plain seconds per step:

0.017057375 0.017202624999999999 0.017105124999999999
0.017302125000000002 0.017093625000000001 0.017236749999999999

Live seconds per step:

0.17001749999999999 0.16580125000000001 0.16527125000000001
0.16262499999999999 0.16294375 0.16808999999999999

Postprocess wall seconds:

0.50932897700113244 0.52051570800540503 0.51081940299627604
0.50154747300257441 0.50169597800413612 0.50100433899933705

## lammps 64000 waters every 10

Warmup plain seconds per step:

0.017345300000000001

Warmup live seconds per step:

0.099632699999999991

Plain seconds per step:

0.0171755 0.018377600000000001 0.016941400000000002 0.017250700000000001
0.017012599999999999 0.019087200000000002

Live seconds per step:

0.0931588 0.086526999999999993 0.087532799999999994 0.086898100000000006
0.088033100000000003 0.104407

Postprocess wall seconds:

0.62307991999841761 0.52231197400396923 0.52448624800308608
0.51828424699488096 0.52810284900624538 0.5357505259962636

## lammps 64000 waters every 100

Warmup plain seconds per step:

0.023650600000000001

Warmup live seconds per step:

0.032671899999999997

Plain seconds per step:

0.023604300000000002 0.0240748 0.0237418 0.023469000000000004
0.023447100000000002 0.023523299999999997

Live seconds per step:

0.031735300000000001 0.032543900000000001 0.031431000000000001
0.031352600000000001 0.032055300000000002 0.031885200000000002

Postprocess wall seconds:

0.51479788700089557 0.51772812900162535 0.51025588500488084
0.56786025700421305 0.52662060500006191 0.5030123780015856

# GROMACS native module size and cadence

The same diamond lattices as the LAMMPS table, one site per water.
Plain is mdrun.
Live is the d-SEAMS module compiled into mdrun.
Both times are the mdrun wall time divided by the number of steps.
The run uses one rank and one thread.
The analysed-frame column is the median extra seconds per step
multiplied by the cadence.

| Waters | Cadence | Steps | Median absolute | Median ratio | Per analysed frame |
| --- | --- | --- | --- | --- | --- |
| 1000 | 1 | 40 | 0.0013000000000000002 | 9.6666666666666661 | 0.0013000000000000002 |
| 1000 | 10 | 40 | 0.00037500000000000006 | 3.5 | 0.0037500000000000007 |
| 1000 | 100 | 100 | 0.00021000000000000004 | 2.5 | 0.021000000000000005 |
| 64000 | 1 | 8 | 0.1255 | 14.199473684210528 | 0.1255 |
| 64000 | 10 | 10 | 0.051650000000000001 | 5.8078078078078068 | 0.51649999999999996 |
| 64000 | 100 | 100 | 0.012189999999999999 | 2.546969689649452 | 1.2189999999999999 |

# Analysis time against the command line

The in-situ column is the extra loop time of one analysed frame.
The command line included column is the wall time of seams rings,
chill-plus and cages.
Process start and file reading are included.
The excluded column is those three censuses after the dump is already in
memory.

| System | In-situ per frame | CLI included | CLI excluded |
| --- | --- | --- | --- |
| 1000 waters, one frame | 0.0015129174999999999 | 0.026832789997570217 | 0.0086144029999999996 |

The in-situ number is the median extra time per step of the 1000-water
every-step run (40 steps).
It includes the first frame.
The command line excluded time is one frame.
Its included time is one process per census.

# OpenMP threads

The ladder times a warm cached workload.
Five dumped frames are parsed before the clock starts.
One engine stays open for the warmup and the six samples.
Each sample calls observe on those frames and divides by five.
Parsing, construction, and the cold first pass stay outside the samples.
A larger thread count scales when its third quartile is below the
previous first quartile.

Waters 1000, threads 1.
First quartile 0.0011330264998832718.
Median 0.0011367245002475101.
Third quartile 0.001140958000542014.
Waters 1000, threads 2.
First quartile 0.00091119405042263681.
Median 0.00092554080038098618.
Third quartile 0.00093689580025966277.
Waters 1000, threads 4.
First quartile 0.00083847745008824859.
Median 0.00085096440016059201.
Third quartile 0.00085889990041323485.
Waters 1000, threads 8.
First quartile 0.00079719825043866877.
Median 0.00081087250000564379.
Third quartile 0.00081489874901308212.
Waters 64000, threads 1.
First quartile 0.071300232749854323.
Median 0.071487236599932655.
Third quartile 0.072266322999348639.
Waters 64000, threads 2.
First quartile 0.053399939500013711.
Median 0.053623207700729834.
Third quartile 0.053720375249758943.
Waters 64000, threads 4.
First quartile 0.041238538100515143.
Median 0.041474410799855838.
Third quartile 0.041695961650111704.
Waters 64000, threads 8.
First quartile 0.039573009550440477.
Median 0.039789857999858208.
Third quartile 0.039873184550015144.

## lammps raw samples

Plain lammps loop, seconds per step.

```
0.00046736750000000006
0.00045604749999999996
0.00046040999999999996
0.00045836
0.00048961999999999996
0.00045961999999999999
```

In-situ lammps loop, seconds per step.

```
0.0020771775000000001
0.0020220849999999999
0.0019486525000000001
0.0019856624999999998
0.0019509925000000001
0.0019581525000000001
```

Extra lammps time, seconds per step.

```
0.0016098100000000001
0.0015660374999999999
0.0014882425
0.0015273024999999996
0.0014613725
0.0014985325000000001
```

## gromacs raw samples

Plain gromacs loop, seconds per step.

```
0.00015000000000000001
0.00015000000000000001
0.00015000000000000001
0.00015000000000000001
0.00015000000000000001
0.00015000000000000001
```

In-situ gromacs loop, seconds per step.

```
0.001475
0.0014500000000000001
0.0015
0.0014500000000000001
0.0014500000000000001
0.0014500000000000001
```

Extra gromacs time, seconds per step.

```
0.001325
0.0013000000000000002
0.0013500000000000001
0.0013000000000000002
0.0013000000000000002
0.0013000000000000002
```

## Command clock

Six excluded-clock samples, seconds, for one 1000-water frame.

```
0.0099659569999999992
0.0086238600000000006
0.0086286850000000005
0.0086049460000000005
0.0085713319999999992
0.0085551759999999994
```

## Engine match

OpenMM matched 8 frames, including its disk coordinate dump.
One LAMMPS plugin rank reproduced 8 frames.
The two-rank plugin run kept that count.
The in-tree package count was 8.
GROMACS matched 9 frames from its trajectory file.

## LAMMPS routes

One thousand waters, one rank, stride 10, 40 steps.
Six samples follow one warmup.
Plugin seconds per step:

```
0.00090272000000000008
0.00090465249999999999
0.000897185
0.00089108250000000003
0.0010263675
0.00092189750000000003
```

Package seconds per step:

```
0.00088538000000000006
0.00089613750000000004
0.00089916249999999996
0.00088973999999999998
0.00088244750000000002
0.00090023749999999997
```

The plugin first quartile is 0.00089856875000000002.
Its third quartile is 0.00091758624999999994.
The package first quartile is 0.0008864700000000001.
Its third quartile is 0.00089840625000000001.
The package third quartile is below the plugin first quartile.

