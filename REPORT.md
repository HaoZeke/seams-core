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
Measured commit `264f90bb0a70d86b283678da951a983ad5ed6c21`.
Measured tree `5a3b95ba0bac037c619800002d427f3d6040c03c`.
libyodaLib SHA-256
`ca2e511ec114c4118209be7b9d9180bd21424dc6281f13c8d979862523fe1a04`.
LAMMPS SHA-256
`2ec47774f5ecd620929f4076ca033bbffec41a66f07331099ab3d0009f4cc5ad`.
GROMACS SHA-256
`dac429749a5a1dd5eb7c7c9249b59d9fd88cd454f46a959f7b5be3f53690020a`.
Plugin SHA-256
`1f6897cd7b105289877422ee24352d62ac70ea9638c881ed5179c07319ff7f93`.

| Code | Steps | Frames | Q1 overhead | Median | Q3 | IQR | Median ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lammps | 40 | 40 | 0.001554459375 | 0.0015897062500000001 | 0.0016259112500000002 | 7.1451875000000126e-05 | 4.3472568473675199 |
| gromacs | 40 | 40 | 0.0013062500000000001 | 0.0013500000000000001 | 0.0013749999999999999 | 6.8749999999999801e-05 | 9.9999999999999982 |

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
| 1000 | 1 | 40 | 0.0015897062500000001 | 4.3472568473675199 | 0.0015897062500000001 | 0.024714980001590448 | 114451 |
| 1000 | 10 | 40 | 0.00046874124999999998 | 2.0004003880712613 | 0.0046874124999999999 | 0.02366294900275534 | 114451 |
| 1000 | 100 | 100 | 0.00015106400000000004 | 1.3283951058459671 | 0.015106400000000004 | 0.024231931500253268 | 114451 |
| 64000 | 1 | 8 | 0.17781649999999999 | 10.284829497412645 | 0.17781649999999999 | 0.56083858450074331 | 7546706 |
| 64000 | 10 | 10 | 0.071466100000000005 | 5.1115963154166213 | 0.71466099999999999 | 0.53723326599720167 | 7546706 |
| 64000 | 100 | 100 | 0.0083227999999999983 | 1.3423146294227157 | 0.8322799999999998 | 0.53248352450100356 | 7546706 |

## lammps 1000 waters every 1

Warmup plain seconds per step:

0.00047321000000000006

Warmup live seconds per step:

0.0022968800000000003

Plain seconds per step:

0.00046789500000000002 0.00047076749999999998 0.00048377749999999997
0.00053470250000000003 0.00048220000000000001 0.00050893750000000006

Live seconds per step:

0.0020602100000000003 0.0021078774999999999 0.0023939625000000001
0.0020782825 0.0020692975000000001 0.001958245

Postprocess wall seconds:

0.025963620995753445 0.024366796998947393 0.024620978998427745
0.024889469998015556 0.023339667000982445 0.02480898100475315

## lammps 1000 waters every 10

Warmup plain seconds per step:

0.00051006999999999995

Warmup live seconds per step:

0.00090487

Plain seconds per step:

0.00049265500000000005 0.00045432750000000004 0.00045715749999999998
0.00045870249999999997 0.00050096500000000003 0.0004795725

Live seconds per step:

0.00091426000000000007 0.00098978499999999997 0.00090679999999999992
0.00090700999999999998 0.0010148225 0.00096741250000000002

Postprocess wall seconds:

0.026302432001102716 0.023726939005427994 0.022859639997477643
0.023655321005207952 0.023670577000302728 0.023321989996475168

## lammps 1000 waters every 100

Warmup plain seconds per step:

0.00045669300000000004

Warmup live seconds per step:

0.00062414100000000004

Plain seconds per step:

0.00046066999999999999 0.00045652799999999999 0.00045841399999999998
0.00046688699999999998 0.00045468500000000003 0.00047547399999999998

Live seconds per step:

0.00061140300000000005 0.000619201 0.00060950100000000003
0.00061792800000000001 0.00061239699999999998 0.00061349099999999997

Postprocess wall seconds:

0.027290822996292263 0.024707813005079515 0.02255074999993667
0.023707622996880673 0.024159330001566559 0.024304532998939976

## lammps 64000 waters every 1

Warmup plain seconds per step:

0.017572125000000001

Warmup live seconds per step:

0.18073375

Plain seconds per step:

0.017333000000000001 0.018206375 0.025037875000000001
0.017922625000000001 0.019391874999999999 0.020730624999999999

Live seconds per step:

0.17258875000000001 0.20227875000000001 0.18735625 0.1902025 0.2116275
0.20408375000000001

Postprocess wall seconds:

0.56186535800225101 0.57982306399935624 0.55965144400397548
0.5598118109992356 0.59161952899739845 0.52500101199984783

## lammps 64000 waters every 10

Warmup plain seconds per step:

0.018598400000000001

Warmup live seconds per step:

0.091312699999999997

Plain seconds per step:

0.0169832 0.017288500000000002 0.0175704 0.019019000000000001
0.017317199999999998 0.018037900000000003

Live seconds per step:

0.088366399999999998 0.089054099999999997 0.089119400000000001
0.088573100000000002 0.092868800000000001 0.0885986

Postprocess wall seconds:

0.53891432999807876 0.53605979299754836 0.53644634999363916
0.52594016700459179 0.55269166900689015 0.53802018200076418

## lammps 64000 waters every 100

Warmup plain seconds per step:

0.024103900000000001

Warmup live seconds per step:

0.032514599999999998

Plain seconds per step:

0.023865500000000001 0.026931300000000002 0.026262199999999999
0.023781500000000001 0.0247297 0.023909699999999999

Live seconds per step:

0.032562500000000001 0.0323752 0.032630599999999996 0.034615300000000002
0.033062999999999995 0.032222000000000001

Postprocess wall seconds:

0.53274416099884547 0.53897759899700759 0.54993313599698013
0.53222288800316164 0.52956354699563235 0.52572425600374117

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
| 1000 | 1 | 40 | 0.0013500000000000001 | 9.9999999999999982 | 0.0013500000000000001 |
| 1000 | 10 | 40 | 0.00038749999999999999 | 3.583333333333333 | 0.003875 |
| 1000 | 100 | 100 | 0.00023999999999999998 | 2.6547619047619047 | 0.023999999999999997 |
| 64000 | 1 | 8 | 0.13425000000000001 | 14.897299578059073 | 0.13425000000000001 |
| 64000 | 10 | 10 | 0.055700000000000006 | 5.883934169278997 | 0.55700000000000005 |
| 64000 | 100 | 100 | 0.012950000000000001 | 2.6107401518068993 | 1.2950000000000002 |

# Analysis time against the command line

The in-situ column is the extra loop time of one analysed frame.
The command line included column is the wall time of seams rings,
chill-plus and cages.
Process start and file reading are included.
The excluded column is those three censuses after the dump is already in
memory.

| System | In-situ per frame | CLI included | CLI excluded |
| --- | --- | --- | --- |
| 1000 waters, one frame | 0.0015897062500000001 | 0.02470720150085981 | 0.0090290070000000004 |

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
First quartile 0.0011323635506414577.
Median 0.0011494085003505459.
Third quartile 0.0011720731996319957.
Waters 1000, threads 2.
First quartile 0.00092211420087551232.
Median 0.0009318172007624526.
Third quartile 0.00094214810014818789.
Waters 1000, threads 4.
First quartile 0.00081501020067662453.
Median 0.00087310029994114302.
Third quartile 0.00089599304919829592.
Waters 1000, threads 8.
First quartile 0.00078178295007091947.
Median 0.00078593370053567921.
Third quartile 0.0008025637001992437.
Waters 64000, threads 1.
First quartile 0.07167232655083354.
Median 0.07182950509959482.
Third quartile 0.072205959499115124.
Waters 64000, threads 2.
First quartile 0.056023931249728778.
Median 0.056254555899795375.
Third quartile 0.056645876449329086.
Waters 64000, threads 4.
First quartile 0.041969986900221558.
Median 0.043447067800298098.
Third quartile 0.046613623149823979.
Waters 64000, threads 8.
First quartile 0.037015507200339928.
Median 0.037312960200506493.
Third quartile 0.037538295749982353.

## lammps raw samples

Plain lammps loop, seconds per step.

```
0.00046789500000000002
0.00047076749999999998
0.00048377749999999997
0.00053470250000000003
0.00048220000000000001
0.00050893750000000006
```

In-situ lammps loop, seconds per step.

```
0.0020602100000000003
0.0021078774999999999
0.0023939625000000001
0.0020782825
0.0020692975000000001
0.001958245
```

Extra lammps time, seconds per step.

```
0.0015923150000000002
0.00163711
0.001910185
0.0015435800000000001
0.0015870975000000002
0.0014493074999999999
```

## gromacs raw samples

Plain gromacs loop, seconds per step.

```
0.00015000000000000001
0.00015000000000000001
0.00015000000000000001
0.000175
0.00015000000000000001
0.00015000000000000001
```

In-situ gromacs loop, seconds per step.

```
0.0015249999999999999
0.001475
0.0018249999999999998
0.001475
0.0014500000000000001
0.0015249999999999999
```

Extra gromacs time, seconds per step.

```
0.0013749999999999999
0.001325
0.0016749999999999998
0.0012999999999999999
0.0013000000000000002
0.0013749999999999999
```

