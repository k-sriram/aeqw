# Automatic Equivalent Width

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Gitmoji](https://img.shields.io/badge/gitmoji-%20😜%20😍-FFDD67.svg?style=flat-square)](https://gitmoji.dev)

Copyright © 2017 Sriram Krishna

## Introduction

`aeqw` is a python script which can be used to calculate the abundance of an element in a stellar photosphere which would produce a prescribed amount of equivalent width for a particular spectral line. This is used in conjunction with the code `SYNSPEC`. `SYNSPEC` simulates radiative transfer through a given model photosphere. Among other things it calculates the equivalent linewidths for a given spectral line in the absorption spectrum that is expected to emerge from the photosphere. This is designed for use with `SYNSPEC 51` that comes along with `TLUSTY 205`.

This document begins with a [brief overview](#how-to-use-synspec) of `SYNSPEC`. Then explains [how to use `aeqw`](#how-to-use-aeqw). Finally it describes [how `aeqw` functions](#how-aeqw-functions) if one needs to know how to debug/modify the program.

## How to use `SYNSPEC`

For a detailed explanation of how to use `SYNSPEC` read [arXiv:1706.185859](https://arxiv.org/abs/1706.01859) § 5. Input/output is done `FORTRAN` unit files, on unix systems this should be of the form `fort.n`, where n is an integer. We shall discuss only the details `aeqw` concerns with. These are `fort.19`, `fort.55`, and `fort.56` for input and `fort.16` for output. The model information is as generated from the `TLUSTY` program. If the name of the model is `modelname`, then the input parameters will be in `modelname.5` and the model photosphere in `modelname.7`. To run the program use:

```sh
RSynpec modelname
```

### `fort.19`
This file contains information about all the lines that we wish the program to consider. This file contains one spectral line in each line of the file. In the `TLUSTY` `FORTRAN` code this format is named `INLIN`. It looks like:

```
  425.3589 16.02  0.107  147146.000 2.0  170648.938 3.0    0.00   0.00   0.00 0
  428.4979 16.02 -0.233  146736.547 1.0  170067.313 2.0    0.00   0.00   0.00 0
  433.2692 16.02 -0.564  146696.188 0.0  169770.047 1.0    0.00   0.00   0.00 0
  435.4566 16.02 -0.959  147690.984 2.0  170648.938 3.0    0.00   0.00   0.00 0
  436.1527 16.02 -0.606  147146.000 2.0  170067.313 2.0    0.00   0.00   0.00 0
  436.4747 16.02 -0.805  147744.547 3.0  170648.938 3.0    0.00   0.00   0.00 0
  449.9245 16.02 -1.640  147550.313 1.0  169770.047 1.0    0.00   0.00   0.00 0
```

Here’s an explanation of each component (with column numbers and description):

```
[ 0-10] ALAM	- wavelength (in nm)
[11-16] ANUM	- code of the element and ion (as in Kurucz-Peytremann)
		  (eg. 2.00 = HeI; 26.00 = FeI; 26.01 = FeII; 6.03 = C IV)
[17-23] GF	- log gf
[24-35] EXCL	- excitation potential of the lower level (in cm*-1)
[36-39] QL	- the J quantum number of the lower level
[40-51] EXCU	- excitation potential of the upper level (in cm*-1)
[52-55] QU	- the J quantum number of the upper level
[57-63] AGAM	= 0. - radiation damping taken classical
		> 0. - the value of Gamma(rad)
[64-70] GS	= 0. - Stark broadening taken classical
		> 0. - value of log gamma(Stark)
[71-77] GW	= 0. - Van der Waals broadening taken classical
		> 0. - value of log gamma(VdW)
[78-79] INEXT	= 0 - no other record necessary for a given line
		> 0 - next record is read, which contains:
```

For all the lines we use `INEXT = 0`. The data is demarked by columns indices.

### `fort.55`

This files contains various parameter values and switches that are required by `SYNSPEC`. These have to be set according to how you want to generate the required synthetic spectra and a better knowledge of the `SYNSPEC` program would let you know how to use this. Read the `SYNSPEC` documentation to know more. Here's a brief template of the file:

```
imode   idstd   iprin
inmod   intrpl  ichang  ichemc
iophli  nunalp  nunbet  nungam  nunbal
ifreq   inlte   icontl  inlist  ifhe2
ihydpr  ihe1pr  ihe2pr
alam0   alast   cutof0  cutofs  relop   space
nmlist,  (iunitm(i), i=1, nmlist)
vtb
nmu0    ang0    iflux
```

See [arXiv:1706.01859](https://arxiv.org/abs/1706.01859) § 5.3.1 to know what each parameter does in detail. Here is an example `fort.55` file.

```
1 32 0
1 0 0 1
0 0 0 0 0
1 0 0 0 0
0 2 0
4147.5 4157.5 15 50 1.0-09 0.010000
0 0i
40.000000
```

`aeqw` only modifies `alam0`, `alast`, and `relop`. Also `ichemc` should be set to 1 for `aeqw` to function and `space` should be small enough that bins are of appropriate size. The last parameter is the maximum spacing between wavelength points in the spectrum in `fort.7` output. We read the output from `fort.16`. Here various wavelengths are binned together. It appears that approximately this value represents the bin size in units of 10&nbsp;nm.

### `fort.56`

This file is used to set the abundances for elements for which we want the abundance to be different from what is given in the model file. The first line of this file contains the number of elements which would be entered. The following lines contain information about each element. We need to specify the atomic number and abundance. This can be set in three ways:
1. As abundance relative to another atom. Set it as the positive number with value equal to this ratio. The reference element is the one with atomic number `ITLAS` in TLUSTY. By default this is H. This can be changed by changing the `ITLAS` in the non-standard parameter file (filename in thrid line of `modelname.5`)
2. As the solar abundace, for this set it to 0.
3. As a multiple of the solar abundance, for this set it to the negative multiple of the desired value.

Example:
```
1
14 3.45-05
```

If no information is provided it takes the default value from the model file. To test for zero abundance, we use a very small value like 10-10.

### `fort.16`

This output file tells the equivalent linewidth in for each wavelength. The program divides the wavelength space into a large number of bins (the size of these bins is specified in `fort.55`). The output has information about each bin on one line. The first two numbers are the bounds of the bin and the third number is the equivalent width (in mÅ). All the bins near a particular spectral line can be summated to calculate the synthetic equivalent linewidth of the line. Here’s a sample file:

```
    4494.200    4495.367         0.0         0.0         0.0         0.0
    4495.367    4496.533         0.0         0.0         0.0         0.0
    4496.533    4497.700         0.0         0.0         0.0         0.0
    4497.700    4498.866         0.1         0.1         0.1         0.1
    4498.866    4500.034        15.9        15.9        16.0        16.0
    4500.034    4501.203         0.0         0.0         0.0         0.0
    4501.203    4502.372         0.0         0.0         0.0         0.0
    4502.372    4503.541         0.0         0.0         0.0         0.0
    4503.541    4504.200         0.0         0.0         0.0         0.0
```

## How to use `aeqw`

### Installation and running
Have the files `autoeqw.py`, `isynspec.py`, and `aeqw.conf` in desired directory or add a link to `autoeqw.py`. `autoeqw.py` should be executable. If it is not use the following command to set it:

```sh
chmod u+x autoeqw.py
```

To run the program simply type the following in the shell (if the name of the model is `hhe35lt`):

```sh
./autoeqw.py hhe35lt
```

### Before Running

Set up `fort.55` as per your liking. Also have the model file in the directory. Edit the config file (`aeqw.conf`) and input file (default `aeqw.in`) as required.

### Configuration file `aeqw.conf`

The configuration information is stored in `aeqw.conf` in a format similar to the `.ini` format. Copy the example `aeqw.conf` file and modify to your needs. Here is an example configuration file:

```ini
[DEFAULT]
INFN = aeqw.in
OUTFN = aeqw.out
EXTRALOGFN = 
INITABUN = 1e-4
NULLABUN = 1e-10
LOGATREF = 11.54
BROAD = 2.0
RANGE = 5.0
EPSILON = 0.1
SEP19 = False

[aeqw]
# Put your custom configuration here
# Such as 
OUTFN = ciiilines.aeqwout
LOGATREF = 10.50
```
The meaning of these parameters are:  
**`INFN`**: The name of the input file.  
**`OUTFN`**: The name of the output file.  
**`EXTRALOGFN`**: If this is specified the log data is also written to this file apart from `aeqw.log`.  
**`INITABUN`**: The initial assumed abundance in the absence of anything given in `fort.56` while estimating the abundance. See [How `aeqw` functions](#how-aeqw-functions) point 5.iii.  
**`NULLABUN`**: The abundance used to estimate the zero values of equivalent width.  See [How `aeqw` functions](#how-aeqw-functions) point 5.ii.  
**`LOGATREF`**: Logarithm of the absolute abundance of the reference element. This is used to calculate the absolute abundance of all the elements in the output file.  
**`BROAD`**: Half width (in Å) upto which absorption is assumed to come from the line. Set it so as to cover the entire line. If set correctly 'wing%' should be low for isolated lines and non-isolated lines shouldn't feed into each other.  
**`RANGE`**: Half width (in Å) of the generated synthetic spectrum used for analysis. This can be much larger than the linewidth, smaller values just save compute time.  
**`EPSILON`**: Accuracy to which the program will try to match the equivalent width.  

The configuration parameters can be overriden by passing them as command-line arguments. Run the following code to see how to do it.

```sh
./autoeqw.py -h
```

### Input file: `aeqw.in`

An example input file looks like:

```
C
C C iii Ei 323(1)k
  405.6061  6.02  0.267  324212.490 2.0  348859.990 3.0   10.23   0.00   0.00 0 159.4
  415.2514  6.02 -0.112  323076.880 1.0  347151.890 2.0    8.80   0.00   0.00 0 174.0
  415.6504  6.02  0.059  323101.360 2.0  347153.260 3.0    8.80   0.00   0.00 0 199.2
  416.2877  6.02  0.218  323140.330 3.0  347155.410 4.0    8.80   0.00   0.00 0 229.2
  418.6900  6.02  0.918  322702.020 3.0  346579.310 4.0    0.00   0.00   0.00 0 320.0
  524.9112  6.02 -0.316  324212.490 2.0  343258.030 1.0    9.68   0.00   0.00 0  84.3
C
C C iii Ei 309(1)k
  432.5561  6.02 -0.759  310006.320 1.0  333118.210 2.0    8.96   0.00   0.00 0 131.5
  465.9058  6.02 -0.654  308248.910 1.0  329706.470 1.0    9.62   0.00   0.00 0 145.7
  466.3642  6.02 -0.530  308248.910 1.0  329685.380 0.0    9.62   0.00   0.00 0 114.9
  466.5860  6.02  0.044  308317.290 2.0  329743.570 2.0    9.62   0.00   0.00 0 273.4
  467.3953  6.02 -0.433  308317.290 2.0  329706.470 1.0    9.62   0.00   0.00 0 202.6
  524.4665  6.02 -1.183  308216.580 0.0  327278.270 1.0    9.71   0.00   0.00 0  56.9
  525.3575  6.02 -0.707  308248.910 1.0  327278.270 1.0    9.72   0.00   0.00 0 121.3
  527.2522  6.02 -0.486  308317.290 2.0  327278.270 1.0    9.72   0.00   0.00 0 163.2

```

Each line specifies information about the spectral lines and their target equivalent width. Blank lines and lines which start with a `#` or a `C` are ignored. Whatever follows a `C` in a line which starts with a `C` is also copied as is to the output file.
The first part of each line uses the same format in `fort.19`. Following this information one can input the target equivalent width in mÅ. This can be entered as a zero in order to combine this line with the following line. This will force the program to calculate the equivalent width for both lines together with the value specified in the next line to be the combined equivalent width. If this number is left blank then nothing is evaluated for this line but is still added to the `fort.19` file.

The estimated initial values of the abundances can be entered in the `fort.56` file in the following format:

```
2
6 5.30-07
7 3.25-06
```

The first line tells that there are two elements, the following lines tell that the abundance of carbon in 5.30×10⁻⁷ and that of nitrogen is 3.25×10⁻⁶.

### Output file: `aeqw.out`
Here is the output file that the example input file will generate.

```
hhe35lt 33000.00 4.00
LAMBDANM   Z.Q   ABUN/ref  LOGABUN   wing%

C iii Ei 323(1)k
405.6061   6.02  1.56e-03     8.73      1%
415.2514   6.02  4.40e-03     9.18    117%
415.6504   6.02  3.93e-03     9.13     82%
416.2877   6.02  3.74e-03     9.11      2%
418.6900   6.02  1.66e-03     8.76      2%
524.9112   6.02  4.33e-03     9.18    109%

C iii Ei 309(1)k
432.5561   6.02  7.30e-03     9.40      0%
465.9058   6.02  7.53e-03     9.42     93%
466.3642   6.02  1.34e-03     8.67    100%
466.5860   6.02  3.25e-03     9.05     18%
467.3953   6.02  8.68e-03     9.48      1%
524.4665   6.02  8.58e-03     9.47    187%
525.3575   6.02  8.60e-03     9.47     84%
527.2522   6.02  8.84e-03     9.49      1%
```

The first line contains the temperature and logarithm of surface gravity. The next line is just a header line. The blank lines and the lines which read like “Silicon II lines” are generated due to the lines in the input file which start with a `C`. The comments have been added as is to the output file. The other lines contain the output. Each line has the wavelength, the atomic number, the level of ionization -1, the relative abundance, and the logarithm of the absolute abundance.

## How `aeqw` functions

`aeqw` calculates abundances by a method of trial and error (you could call it estimated guesses). It simply writes a guess abundance in `fort.56`, runs the `SYNSPEC` program and then reads the equivalent width from `fort.16`. Then it refines its guess and reiterates till it gets a value of equivalent width close enough to the target value. It also automates the job for a large number of lines. The following lines explain how it works.

1. All the parameters from `fort.55` are read. This is so that the program is able to write to it later.
2. The input file aeqw.in is read. The first line contains a few parameters. The following lines contain information about all the spectral lines to be tested and the target equivalent widths. Some lines which are too close to each other can be grouped and a common equivalent width be given (See the specifications in the next section to know how to enter input in this file). The data is stores as an array of sets of lines. Each set contains lines which need to be evaluated together. Most sets will contain single lines.
3. It is checked if `fort.8` and the model file have the same temperature and surface gravity. If not a warning is displayed. However, the execution is not interrupted.
4. All lines are written to `fort.19`.
5. The following steps are taken for each set of lines:
   1. The bounds of the synthetic spectrum are set. This is simply (smallest wavelength – `RANGE`, largest wavelength + `RANGE`) where `RANGE` is a parameter that can be set in `aeqw.conf` (this can be much larger than the linewidth, smaller values just save compute time). This is written to `fort.55`.
   2. The equivalent width is calculated (using the method in steps iv. to vi.) for an abundance of `10e-10` (settable by modifying `NULLABUN` in `aeqw.conf`). This is used as a zero baseline for future calculations.
   3. An initial value of abundance is assumed. (settable by modifying `INITABUN` in `aeqw.conf`, default `1e-4`).
   4. The assumed value of abundance is written into `fort.56`. The atomic number is inferred from the specification of the line.
   5. `SYNSPEC` is run.
   6. `fort.16` is read. The equivalent width is calculated. Only bins which are up to a distance specified by the parameter `BROAD` (set it so that it covers all absorbtion, but not large enough to read from other lines) from the spectral line are considered. If the edge of the considered range is inside a bin, the bin is considered partially.
   7. It is checked if the value of equivalent width is acceptable (using the parameter `EPSILON` in `aeqw.conf`). If not a new estimate for abundance is made and the steps iv. to vi. are repeated. The new estimate is arrived by assuming the equivalent width to be a linear function of abundance. It is also checked if the line is too weak or if we see emission.
6.  The output is written to the output file. See specifications in [previous section](#how-to-use-aeqw) to interpret it.

## Logging

The program logs each step of the working in the file `aeqw.log`. The last 10 logs are also stores in `aeqw.log.n`.
