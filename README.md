# mri-homomorphic-noise-estimation
Python implementation of “Spatially variant noise estimation in MRI A homomorphic approach”

## How to run
### Requirements
Program was developed using Anaconda python distribution which includes regular python interpreter and additional popular packages for sience, math, data analysis etc.

Anaconda can be downloaded from http://continuum.io/downloads. After installing environment is ready to be used and run application.

In thoery program should run on any 2.7.x python distribution with `scipy` and `numpy` packages installed but it wasn't tested in that way.

### Command line
Program can be executed from command line in a following way. Make sure `python` command points to distribution from Anaconda.

```
$ python run.py -e            # runs algorithm with predefined sample input
$ python run.py -c config.ini # runs algorithm with parameters defined in config.ini file
```

Options `-e` and `-c` are mutually exclusive.

### Output
Program produce noise maps in csv format in `/data/output` directory.

## Versions
Software:
* Version of Anaconda used: 2.2.0 (64-bit)
* Version of Python used: 2.7.9 (64-bit)

Tested on:
* Windows 7 Home Premium with Service Pack 1 (64-bit)

## Known issues
None
