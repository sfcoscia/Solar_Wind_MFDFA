# PSP & MMS1 Solar Wind MFDFA 
This repo contains the MFDFAlib.py which acts as a library of Python code that
performs Multifractal De-trended Fluctuation Analysis on MMS1
and PSP magnetic field time series data at sub-Alfvenic and super-Alfvenic
intervals for comparison. The analyses are performed in Jupyter notebooks. 

### Installation

First, clone the repo onto your local machine. 

#### Prerequisites 
Python 3 or later must be installed. One also needs a way of viewing a Jupyter
Notebook.

Install the required libraries:
```python
pip3 install requirements.txt
``` 

### MMS1 Intervals
To access the MMS1 interval analysis, enter the MMS1 Interval directory:
```bash
cd MMS1_Intervals
```
Access the Jupyter notebook "MFDFA_MMS1.ipynb" with the environment of choice for the MFDFA plots and code. For
example, with VS Code run:
```bash
code MFDFA_MMS1.ipynb
```

Example plots are generated in the "plots" directory. 

### PSP Intervals
From the main project directory, to access the MFDFA of select PSP intervals, 
enter the "PSP_Intervals" directory: 
```bash
cd PSP_Intervals
```
Once in this directory, there are several MFDFA notebooks performed on different
sets of sub and super-Alfvenic intervals from PSP. They are each indexed by
1,2,3,... For example, the first set has the file name
"MFDFA_PSP_Interval1.ipynb".
(Note to self : Change this to the date of the interval)

The "ACF" notebook is a notebook containing the a test of the autocorrelation
function of the 1st PSP interval, which is now integrated ionto each MFDFA
notebook. The "plots" folder contains plots for each set of intervals. 

### What is MFDFAlib.py?
This set of functions essentially acts as an intermediate between
pyspedas for loading data and the MFDFA library, while also containing code to
calculate the multifractal spectrum from the MFDFA fluctuation funcion. 

### Fractional Gaussian Motion Tests
There are two notebooks for in the main project directory for the purpose of
testing the MFDFAlib code with fractional Gaussian motion at a set Hurst
parameter. These are titled "FGM_test.ipynb" and "FGM_test.ipynb". (Note to
self: change these names)
