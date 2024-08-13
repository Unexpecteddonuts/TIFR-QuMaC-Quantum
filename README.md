# Documentation v3

## Imported Libraries
```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import cauchy
from scipy.signal import find_peaks
import os
```

This python script starts off by importing the following libraries:

1. `matplotlib.pyplot`: Used for creating static, animated, and interactive visualizations in Python.
2. `numpy`: A library for the Python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays.
3. `curve_fit`: A function from the SciPy library. It is used to fit some data to a function.
4. `cauchy`: A continuous probability distribution. It is used in many areas of statistics, physical sciences, and engineering.
5. `find_peaks`: A function from scipy library used to find peaks of vectors or matrices.
6. `os`: Provides a way of using operating system dependent functionality to read or write from files, manipulate paths, etc.

## Function Definitions

### `cauchy_fit() function`
```python
def cauchy_fit(x, x0, gamma, A, y0):
    return y0 + A * cauchy.pdf(x - x0, scale=gamma)
```
This function calculates the Cauchy distribution function for a given value of x, and parameters x0, gamma, A, and y0, where `x0` is the location parameter (peak of the distribution), `gamma` is the scale parameter (FWHM of the distribution - Full Width at Half Maximum), `A` is the amplitude of the distribution and `y0` is the baseline of the Cauchy distribution.

### `ext_BW() function`
```python
def ext_BW(freq, ydata):
    ...   
    return properties
```
It obtains sharp frequency bands for a given set of frequencies and corresponding harmonic power data. After identifying the peaks of power, it estimates a set of initial parameters, uses them for curve fitting to Cauchy Distribution, then calculates properties and uncertainties.

## Process Flow

The script first defines a range of power values. It then initializes a list of storage variables. Then, for each power level in the defined range, the script loads the corresponding data file, reads in the frequency and y data, and calculates the properties of a resonator based on this data.

These properties are then appended to the previously defined storage lists. In this loop, the lines `... filepath = f'F:\\Expt Data\\date\\folder\\data_{p}.0_dBm.txt'` and `... data = np.transpose(np.loadtxt(filepath))` are particularly important for each file in filepath, load text data and transform it (transposition).

These calculated properties and fit parameters are then saved to varaibles `fit_param_list` and `Q_with_pow_list` respectively and also appended to the individual storage lists - `kint_list`, `kext_list`, `Qint_list`.

Finally, these data lists are written to text files using numpy's `savetxt()` function with appropriate headers. The outputs would be saved in two ".txt" files named "Lorentian_fit_params.txt" and "Q_factors.txt", located in current working diretory. You can find these files by calling `os.getcwd()` in your python console or script, which will return you current working directory..

## Code structure and explantion

```python
# Defining power range
pmin, pmax = -100, 10
del_p = 10
p_vals = np.arange(pmin, pmax + del_p - 1, del_p)
```
Define a range of power levels from `pmin` to `pmax`, with steps of `del_p` units. A list of power values is generated using `np.arange()`.

```python
# storage lists
kint_list, kext_list, Qint_list, fit_param_list, Q_with_pow_list = [], [], [], [], []
```
Here, five lists are initiated as data storage varibles.

```python
# Looping over power values and filepaths
for p in p_vals:
    filepath = f'F:\\Expt Data\\date\\folder\\data_{p}.0_dBm.txt'
    data = np.transpose(np.loadtxt(filepath))
    freq = data[0]
    ydata = data[1]
```
Loop over each power level, and for each level, build a filepath to the corresponding data file, load the data file and then transpose. The frequency is the first element (`data[0]`) of `data`, ydata is the second element(`data[1]`)

```python
properties = ext_BW(freq, ydata)
```
Call `ext_BW()` function on this data to calculate fit parameters and properties.

Next few lines are about generating lists for further analysing.
```python
fit_param_list.append([...])
Q_with_pow_list.append([...])
kint_list.append(properties['internal_bandwidth'])
kext_list.append(properties['external_bandwidth'])
Qint_list.append(properties['internal_quality_factor'])
```

```python
# Saving the data to text files
np.savetxt(os.path.join(os.getcwd(), "Lorentian_fit_params.txt"), fit_param_list,
           header="Power Freq Bandwidth Internal_BW External_BW Internal_Qfactor")
np.savetxt(os.path.join(os.getcwd(), "Q_factors.txt"), Q_with_pow_list,
           header="Power Freq Internal_Qfactor Internal_QFactor_err External_QFactor External_QFactor_err ")
```
Finally, save the calculated data to two separate text file named "Lorentian_fit_params.txt" and "Q_factors.txt". Headers for the columns of data are provided in string format in the `header` argument of `numpy.savetxt()`.

NOTE:
1. `np.transpose()` transforms rows into columns, and vice versa, for a 2D matrix.
2. `os.getcwd()` returns the path of current working directory.
3. `np.savetxt()` writes a 1D or 2D numpy array to text file, use 'header' argument to specify column names.
4. `os.path.join()` combines one or more path names into a single path.
5. `numpy.arange()` generates an array of evenly spaced values in a specified range.
6. `list.append()` adds an element to the end of a list.
