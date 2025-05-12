# Raman Scattering of Rhenium for Secondary Pressure Calibration

## About
### Single Peak Fit 
This is a Python software that iteratively peak fits Raman (or XRD data) with three different background types  
### Diamond Edge Fit
This is a Python software that will fit diamond edge data using flat, bevel, and tDAC anvils.

## Running The Software &nbsp;![Jupyter Notebook](https://img.shields.io/badge/Jupyter-Notebook-orange?logo=jupyter) ![Python](https://img.shields.io/badge/Python-3.x-blue?logo=python)
This software utilized Jupyter Notebooks, allowing the user to execute individual sections of code interactively. Additionally, the Single Peak Fit is available as a standalone Python script.

## Python Dependencies  
This software requires the following Python libraries:
* Core Libraries:
  * `math` - Mathematical functions
  * `os` - Operating system interaction
  * `csv` - CSV file handeling
  * `datetime` - Date and time utilities
* Scientific Computing and Data Processing:
  * `numpy` - Numerical computations (`loadtxt` for loading data)
  * `scipy` - Scientific computing (`curve_fit` for curve fitting)
  * `wofz` - Voigt function from `scipy.special`
  * `derivative` - Numerical differentiation
  * `savgol_filter` - Savizky-Folay filter for smoothing
* Data Visualization:
  * `matplotlib` - Plotting library (`matplotlib.pyplot` for visualization)

## Authors
Allison M. Pease<sup>1, 2</sup>, Claire Zurkowski<sup>1</sup>, Stella Chariton<sup>3</sup>, Heidi N. Krauss<sup>2</sup>, Daniel Sneed<sup>1</sup>, Vitali Prakapenka<sup>3</sup>, Bruce Baer<sup>1</sup>, Susannah M. Dorfman<sup>2</sup>, and Earl F. O'Bannon III<sup>1</sup>  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<sup>1</sup>Physics Division, Physical and Life Sciences Directorate, Lawrence Livermore National Laboratory, Livermore, California, United States  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<sup>2</sup>Department of Earth and Environmental Sciences, Michigan State University, East Lansing, Michigan, United States  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<sup>3</sup>Center for Advanced Radiation Sources, The University of Chicago, Chicago, IL, United States  

## Contact
For questions regarding this software, contact:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Allison M. Pease: peaseall@msu.edu  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Heidi N. Krauss: Heidi.N.Krauss@gmail.com

## License &nbsp;[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
This project is licensed under the GNU General Public License v3.0. See the LICENSE file for details.  

##
Copyright &copy; 2025 Pease et al.

Last Updated: March 25<sup>th</sup>, 2025
