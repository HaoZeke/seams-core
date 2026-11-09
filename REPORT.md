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
Measured commit `f3c648bf21a2c9fad084702f48605714e6e2881e`.
Measured tree `548a9ca3284a98e53949ac9f72802991ecbd1845`.
libyodaLib SHA-256
`de5b160bc827eceb53aea7c7393da5803d2517a0fd21c0e4abe316b4d32aa8c4`.
LAMMPS SHA-256
`2ec47774f5ecd620929f4076ca033bbffec41a66f07331099ab3d0009f4cc5ad`.
GROMACS SHA-256
`6c912b0dc3883370722066345895ee6a128a1d741db3b64c3454408088b0318c`.
Plugin SHA-256
`1f6897cd7b105289877422ee24352d62ac70ea9638c881ed5179c07319ff7f93`.

| Code | Steps | Frames | Q1 overhead | Median | Q3 | IQR | Median ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lammps | 40 | 40 | 0.00150836875 | 0.0015276649999999999 | 0.0016411574999999998 | 0.00013278874999999983 | 4.3750375189641275 |
| gromacs | 40 | 40 | 0.0066812499999999997 | 0.0067999999999999996 | 0.0068249999999999995 | 0.00014374999999999978 | 12.599999999999998 |

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
| 1000 | 1 | 40 | 0.0015276649999999999 | 4.3750375189641275 | 0.0015276649999999999 | 0.024076607500319369 | 114451 |
| 1000 | 10 | 40 | 0.0043239287499999991 | 10.410809417054672 | 0.043239287499999987 | 0.094803200001479127 | 114451 |
| 1000 | 100 | 100 | 0.0027208535000000002 | 6.9910162595030947 | 0.27208535 | 0.071969930995692266 | 114451 |
| 64000 | 1 | 8 | 0.47819706249999999 | 28.882504336957169 | 0.47819706249999999 | 1.5177708290029841 | 7546706 |
| 64000 | 10 | 10 | 0.067820450000000004 | 5.0324964911493471 | 0.6782045000000001 | 0.50354381500073941 | 7546706 |
| 64000 | 100 | 100 | 0.0080616500000000018 | 1.344039522959986 | 0.80616500000000013 | 0.50741066550472169 | 7546706 |

## lammps 1000 waters every 1

Warmup plain seconds per step:

0.00045361749999999997

Warmup live seconds per step:

0.0019494325

Plain seconds per step:

0.00045323000000000005 0.00045601749999999998 0.00045369249999999996
0.00045205499999999998 0.00044991500000000006 0.000450765

Live seconds per step:

0.001967145 0.0021493375 0.0019602125 0.0019934699999999998
0.0021243199999999999 0.0019389175000000002

Postprocess wall seconds:

0.026318415999412537 0.023054285004036501 0.022402334005164448
0.024185517999285366 0.024358021997613832 0.023967697001353372

## lammps 1000 waters every 10

Warmup plain seconds per step:

0.00045348999999999997

Warmup live seconds per step:

0.0049337749999999996

Plain seconds per step:

0.00046347249999999999 0.00045269750000000001 0.00045562000000000001
0.00045722000000000004 0.00046090250000000002 0.00045917500000000003

Live seconds per step:

0.0047338249999999997 0.0045033 0.0048331249999999997
0.0045530500000000003 0.0049410999999999995 0.0049382250000000001

Postprocess wall seconds:

0.099878818000433967 0.1002049220041954 0.069941167996148579
0.080347408002126031 0.089727582002524287 0.10397819400532171

## lammps 1000 waters every 100

Warmup plain seconds per step:

0.00045343

Warmup live seconds per step:

0.00310188

Plain seconds per step:

0.00045371500000000004 0.00045333499999999997 0.00045482799999999995
0.00045579700000000003 0.00045449999999999999 0.00045252599999999997

Live seconds per step:

0.0031894599999999999 0.0030939299999999999 0.00313635
0.0031766100000000003 0.0033803100000000001 0.0031734200000000001

Postprocess wall seconds:

0.10008488300081808 0.071943411996471696 0.076183116005267948
0.071996449994912837 0.064994901993486565 0.06915524099895265

## lammps 64000 waters every 1

Warmup plain seconds per step:

0.017156874999999999

Warmup live seconds per step:

0.67796875000000001

Plain seconds per step:

0.0172795 0.017382499999999999 0.01764775 0.017215999999999999
0.017084874999999999 0.017356125

Live seconds per step:

0.67654124999999998 0.57550000000000001 0.49172874999999999
0.49721874999999999 0.49347625000000001 0.49302499999999999

Postprocess wall seconds:

1.5271113409980899 1.5026323330021114 1.5261975200046436
1.5093441380013246 1.496989028993994 1.5303928789944621

## lammps 64000 waters every 10

Warmup plain seconds per step:

0.016879700000000001

Warmup live seconds per step:

0.085730300000000009

Plain seconds per step:

0.016780300000000001 0.016870099999999999 0.016859800000000001
0.016748599999999999 0.016821900000000001 0.016868399999999999

Live seconds per step:

0.084404300000000002 0.085475099999999998 0.084296599999999999
0.084292100000000009 0.084849800000000003 0.084885299999999997

Postprocess wall seconds:

0.51036503299837932 0.50418833600269863 0.50733081700309413
0.50289929399878019 0.49824416999763343 0.49510163500235649

## lammps 64000 waters every 100

Warmup plain seconds per step:

0.023513700000000002

Warmup live seconds per step:

0.031730299999999996

Plain seconds per step:

0.023440599999999999 0.023423799999999998 0.0238908 0.0238213
0.023500899999999998 0.023536800000000004

Live seconds per step:

0.031636499999999998 0.031351200000000003 0.031375899999999998 0.0316513
0.041307999999999997 0.0319124

Postprocess wall seconds:

0.52784393100591842 0.52164284999889787 0.50448621899704449
0.50604589400609257 0.50283697900158586 0.50877543700335082

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
| 1000 | 1 | 40 | 0.0067999999999999996 | 12.599999999999998 | 0.0067999999999999996 |
| 1000 | 10 | 40 | 0.002 | 4.8095238095238093 | 0.02 |
| 1000 | 100 | 100 | 0.00094499999999999988 | 2.857554994590696 | 0.094499999999999987 |
| 64000 | 1 | 8 | 0.24681250000000002 | 17.184426229508198 | 0.24681250000000002 |
| 64000 | 10 | 10 | 0.052600000000000001 | 6.7649904443382702 | 0.52600000000000002 |
| 64000 | 100 | 100 | 0.01227 | 2.8932125612546877 | 1.2269999999999999 |

# Analysis time against the command line

The in-situ column is the extra loop time of one analysed frame.
The command line included column is the wall time of seams rings,
chill-plus and cages.
Process start and file reading are included.
The excluded column is those three censuses after the dump is already in
memory.

| System | In-situ per frame | CLI included | CLI excluded |
| --- | --- | --- | --- |
| 1000 waters, one frame | 0.0015276649999999999 | 0.023787507499946514 | 0.0086447290000000003 |

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
First quartile 0.0012102289489121175.
Median 0.0012213346002681647.
Third quartile 0.0012299597005039685.
Waters 1000, threads 2.
First quartile 0.00090899695023836111.
Median 0.00091581110027618708.
Third quartile 0.00092321550000633583.
Waters 1000, threads 4.
First quartile 0.00078689519905310594.
Median 0.00079219929975806738.
Third quartile 0.00080019515044114095.
Waters 1000, threads 8.
First quartile 0.0007246010503877188.
Median 0.00072851270015235055.
Third quartile 0.00073041614959947778.
Waters 64000, threads 1.
First quartile 0.071039519449914212.
Median 0.071477220000087993.
Third quartile 0.073258544200507467.
Waters 64000, threads 2.
First quartile 0.053325099450012209.
Median 0.053509692300576714.
Third quartile 0.053626567500032252.
Waters 64000, threads 4.
First quartile 0.041316414350148989.
Median 0.041959418800252023.
Third quartile 0.04281772124959389.
Waters 64000, threads 8.
First quartile 0.036103209249995417.
Median 0.036503266799991246.
Third quartile 0.036542571651079922.

## lammps raw samples

Plain lammps loop, seconds per step.

```
0.00045323000000000005
0.00045601749999999998
0.00045369249999999996
0.00045205499999999998
0.00044991500000000006
0.000450765
```

In-situ lammps loop, seconds per step.

```
0.001967145
0.0021493375
0.0019602125
0.0019934699999999998
0.0021243199999999999
0.0019389175000000002
```

Extra lammps time, seconds per step.

```
0.001513915
0.0016933200000000001
0.0015065199999999999
0.0015414149999999998
0.0016744049999999999
0.0014881525000000001
```

## gromacs raw samples

Plain gromacs loop, seconds per step.

```
0.000125
0.00072500000000000006
0.00052500000000000008
0.00072500000000000006
0.00052500000000000008
0.00052500000000000008
```

In-situ gromacs loop, seconds per step.

```
0.0014
0.0074999999999999997
0.0074749999999999999
0.0073749999999999996
0.0073499999999999998
0.0073499999999999998
```

Extra gromacs time, seconds per step.

```
0.0012750000000000001
0.0067749999999999998
0.0069499999999999996
0.0066499999999999997
0.0068249999999999995
0.0068249999999999995
```


# Build pin

The GROMACS module tip is bdf53dcf3ca265a62600e1e3d5345476d666c6d0.
Upstream main at the branch point is
1ede0115a1c6cc063813d2b13085a86640e4273e.
The module tip is five commits after that branch point.
The clean binary prints 2027.0-dev-20261108-bdf53dcf3c.
The pin is the module tip.
The binary records that commit.

# Frames

OpenMM writes a frame every 5 steps.
The first saved step is 5.
The steps are 5, 10, 15, 20, 25, 30, 35 and 40.
That set has 8 frames.
GROMACS writes a frame at step 0 and then every 5 steps.
The last step is 40.
The steps are 0, 5, 10, 15, 20, 25, 30, 35 and 40.
That set has 9 frames.
The live step list equals the native step list.
OpenMM live and the DCD both use those 8 steps.
GROMACS live and the coordinate file both use those 9 steps.

# Timing head

The size table and the clock were measured at
f3c648bf21a2c9fad084702f48605714e6e2881e.
No file under src changed between that commit and the tip.
The library SHA-256 is the same at both commits.
The live driver then split the ladder helper.
The ladder above was measured again on that helper.
A larger thread count passes when its third quartile is below the
previous first quartile.
The 1000-water ladder and the 64000-water ladder both pass.
