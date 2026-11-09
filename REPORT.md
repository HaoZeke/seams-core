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
Measured commit `8c4657ca26c58e35aa089513c7caccc9100f3779`.
Measured tree `3db38834f3e23b29574edb319b772626bc453d9a`.
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
| lammps | 40 | 40 | 0.001445324375 | 0.0014870025 | 0.00149306875 | 4.7744374999999974e-05 | 4.2263252155432873 |
| gromacs | 40 | 40 | 0.0012812500000000001 | 0.0013000000000000002 | 0.0013000000000000002 | 1.8750000000000103e-05 | 9.6666666666666661 |

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
| 1000 | 1 | 40 | 0.0014870025 | 4.2263252155432873 | 0.0014870025 | 0.024964869498944608 | 114451 |
| 1000 | 10 | 40 | 0.00044709375000000001 | 1.9889480952457042 | 0.0044709375000000001 | 0.026790096002514474 | 114451 |
| 1000 | 100 | 100 | 0.00015411950000000001 | 1.3370893174140113 | 0.015411950000000001 | 0.024794995002594078 | 114451 |
| 64000 | 1 | 8 | 0.148049875 | 9.6133865863196899 | 0.148049875 | 0.51051799000197207 | 7546706 |
| 64000 | 10 | 10 | 0.068839399999999995 | 5.0168142806763516 | 0.68839399999999995 | 0.51757203300076071 | 7546706 |
| 64000 | 100 | 100 | 0.0081123500000000025 | 1.3462169635029637 | 0.81123500000000026 | 0.51200928899925202 | 7546706 |

## lammps 1000 waters every 1

Warmup plain seconds per step:

0.00050197749999999998

Warmup live seconds per step:

0.0019426974999999999

Plain seconds per step:

0.00046272750000000002 0.00045543249999999998 0.00045988000000000002
0.00063904249999999999 0.00045395249999999999 0.00057273250000000001

Live seconds per step:

0.0019455875000000001 0.0019704925000000001 0.0019535899999999998
0.001949515 0.0019450975 0.0020055450000000001

Postprocess wall seconds:

0.039837331998569425 0.024976280998089351 0.026276911004970316
0.024953457999799866 0.024526378001610283 0.024428047996480018

## lammps 1000 waters every 10

Warmup plain seconds per step:

0.00045364250000000004

Warmup live seconds per step:

0.0010218375

Plain seconds per step:

0.00044966 0.000451385 0.00045542750000000001 0.00045050000000000005
0.00045279750000000002 0.0004587975

Live seconds per step:

0.00091065000000000004 0.00089844000000000002 0.0008998124999999999
0.00090247250000000004 0.00089992999999999996 0.00089771750000000009

Postprocess wall seconds:

0.026709539000876248 0.027778447001765016 0.026947218997520395
0.0244118080008775 0.024700916001165751 0.0268706530041527

## lammps 1000 waters every 100

Warmup plain seconds per step:

0.00046238300000000002

Warmup live seconds per step:

0.00060390599999999997

Plain seconds per step:

0.00045713999999999998 0.00045237100000000002 0.00045727200000000005
0.000453657 0.00045843099999999998 0.00045482500000000001

Live seconds per step:

0.00060956400000000005 0.00060884800000000005 0.00061308700000000001
0.00060449199999999999 0.00060261799999999997 0.000628952

Postprocess wall seconds:

0.028975238994462416 0.026026483996247407 0.024023300000408199
0.024032171000726521 0.025557819004461635 0.021891536001930945

## lammps 64000 waters every 1

Warmup plain seconds per step:

0.01757475

Warmup live seconds per step:

0.17011124999999999

Plain seconds per step:

0.017325500000000001 0.0173285 0.017274749999999998 0.017066999999999999
0.017155750000000001 0.017136874999999999

Live seconds per step:

0.17754875000000001 0.16604125 0.16473125 0.16485875 0.16546374999999999
0.16406750000000001

Postprocess wall seconds:

0.521035263001977 0.51043291300447891 0.50675301099545322
0.50840469999820925 0.51060306699946523 0.51114634200348519

## lammps 64000 waters every 10

Warmup plain seconds per step:

0.017593600000000001

Warmup live seconds per step:

0.093568300000000007

Plain seconds per step:

0.017200899999999998 0.0180903 0.0175552 0.017957799999999999
0.017192900000000001 0.0169102

Live seconds per step:

0.089976 0.086515800000000004 0.088295800000000008 0.085573700000000003
0.086033600000000002 0.0857483

Postprocess wall seconds:

0.51125705699814716 0.50773508400015999 0.50847898800566327
0.52497928799857618 0.52451490400562761 0.52388700900337426

## lammps 64000 waters every 100

Warmup plain seconds per step:

0.023848699999999997

Warmup live seconds per step:

0.031905900000000001

Plain seconds per step:

0.023700499999999999 0.024140700000000001 0.023391799999999997
0.023417400000000001 0.0234702 0.024005600000000002

Live seconds per step:

0.033020999999999995 0.031676200000000002 0.031406700000000003
0.033610099999999997 0.03168 0.031601999999999998

Postprocess wall seconds:

0.51586559199495241 0.50622950500110164 0.50428721700154711
0.50815298600355163 0.52640935399540467 0.56716194399632514

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
| 64000 | 1 | 8 | 0.12725 | 14.40647226173542 | 0.12725 |
| 64000 | 10 | 10 | 0.052049999999999999 | 5.8644859813084116 | 0.52049999999999996 |
| 64000 | 100 | 100 | 0.012214999999999998 | 2.5500642549636954 | 1.2214999999999998 |

# Analysis time against the command line

The in-situ column is the extra loop time of one analysed frame.
The command line included column is the wall time of seams rings,
chill-plus and cages.
Process start and file reading are included.
The excluded column is those three censuses after the dump is already in
memory.

| System | In-situ per frame | CLI included | CLI excluded |
| --- | --- | --- | --- |
| 1000 waters, one frame | 0.0014870025 | 0.026611004497681279 | 0.0089189520000000008 |

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
First quartile 0.0011139679998450449.
Median 0.0011187440002686344.
Third quartile 0.0011375511494406965.
Waters 1000, threads 2.
First quartile 0.00094275609990290832.
Median 0.00095209299979615025.
Third quartile 0.00096565794992784491.
Waters 1000, threads 4.
First quartile 0.0008029653999983567.
Median 0.00081115770008182157.
Third quartile 0.00081961370051431004.
Waters 1000, threads 8.
First quartile 0.00071451774929300873.
Median 0.00071575159963686028.
Third quartile 0.00072007304952421698.
Waters 64000, threads 1.
First quartile 0.070575420699969985.
Median 0.070789924199198134.
Third quartile 0.07081062529941845.
Waters 64000, threads 2.
First quartile 0.05265988770042896.
Median 0.052906963600253223.
Third quartile 0.053127067999230351.
Waters 64000, threads 4.
First quartile 0.040522238199264392.
Median 0.04091631929986761.
Third quartile 0.041573661499569418.
Waters 64000, threads 8.
First quartile 0.035112480249517831.
Median 0.035576035000121922.
Third quartile 0.036312606550200147.

## lammps raw samples

Plain lammps loop, seconds per step.

```
0.00046272750000000002
0.00045543249999999998
0.00045988000000000002
0.00063904249999999999
0.00045395249999999999
0.00057273250000000001
```

In-situ lammps loop, seconds per step.

```
0.0019455875000000001
0.0019704925000000001
0.0019535899999999998
0.001949515
0.0019450975
0.0020055450000000001
```

Extra lammps time, seconds per step.

```
0.0014828600000000001
0.0015150600000000001
0.0014937099999999999
0.0013104725
0.0014911449999999999
0.0014328125000000001
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
0.0014500000000000001
0.0014500000000000001
0.00215
0.0014250000000000001
0.0014250000000000001
0.0014500000000000001
```

Extra gromacs time, seconds per step.

```
0.0013000000000000002
0.0013000000000000002
0.002
0.00125
0.0012750000000000001
0.0013000000000000002
```


## Engine match

OpenMM matched 6 frames.
The LAMMPS plugin matched 4 frames.
The LAMMPS package matched 4 frames.
GROMACS matched 5 frames.
Each count is the number of shared steps whose census matches `seams`.

## LAMMPS routes

One thousand waters, one rank, stride 10, 40 steps.
Six samples follow one warmup.
Plugin samples, seconds per step: 0.00090092000000000004,
0.00090812249999999998, 0.00088275750000000003, 0.0009034275,
0.00090755750000000009, 0.00091188999999999994.
Package samples, seconds per step: 0.00089072999999999995,
0.00089272750000000001, 0.00089070499999999999, 0.00089675750000000004,
0.00088873500000000009, 0.00089046249999999991.
The plugin first quartile is 0.000901546875.
Its third quartile is 0.00090798124999999996.
The package first quartile is 0.00089052312499999988.
Its third quartile is 0.00089222812500000005.
The package third quartile is below the plugin first quartile.
