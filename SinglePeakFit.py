# SinglePeakFit.py

    # This software can sort the first column of the original text files.  
        # If the first row of the first column is less than or equal to zero it will be saved as a sorted file in the Data directory  
        # If the first row of the first column is greater than zero it will be saved as a sorted file in the Pressure Calibration directory  
    # The output is a CSV file; the header contains the file name, time processed, average peak positions, standard deviation of peak positions, peak position, full width half max (FWHM), peak fit type, and background fit type  
    # The user can choose which peak fits (Gaussian, Lorentzian, and Voigt) along with which background fits (Constant, Linear, and Second-Order Polynomial) to analyze  
    # The user can choose to analyze positive and/or negative peak fits  
        # To use this for XRD data turn off the negative peak fit  
    # The user can choose for error messages to appear  
    # &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If the standard deviation of the peak position if greater than the user input value, a check standard deviation (CHECK STDV FOR) message will be printed along with the file name, sign, and standard deviation of the peak position  
    # The user can choose to print the anti-stokes average and standard deviation of the peak position for every peak fit  
    # The user can choose to print graphs of the peak fits saved to the Image directory inside the Data directory  
        # A graph for each peak fit and background combination  
        # A comparison graph of all peak fit and background combinations for each text file

# Copyright (C) 2025 Pease et al.

# For questions regarding this software, contact:
    # Allison M. Pease: peaseall@msu.edu
    # Heidi N. Krauss: Heidi.N.Krauss@gmail.com

#######################################################################

# Importing the necessary libraries
import numpy as np
np.set_printoptions(threshold=np.inf)
from numpy import loadtxt
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import wofz
import csv
from datetime import datetime
on = True
off = False

#######################################################################
# User Input Required
#######################################################################

# Folder and File Names and Locations
    # Path to the original directory inputed by the user (Old_Folder) and the new directory (New_Folder)
Path = "C:/Users" 
    # Name of the original directory inputed by the user (Old_Folder)
Old_Folder = "Input_Data"
    # Name of the new directory (New_Folder)
New_Folder = "Output_Data"
    # Name of the CSV file
        # Remember to add .csv to the end of the file name
CSV_File_Name = "Saved_Data.csv"
# Initial Guess
center_guess=180
# Iteration Values
Convergence_Range = 0.05
Max_Iterations = 800
Check_Standard_Deviation_of_Peak_Position = 1.0
# Turn on or off print features
    # When On: this will make a graph for each peak fit and background combination
Make_Images = on
    # When On: this will make a comparison graph of all peak fit and background combinations for each text file
Make_Comparison_Graph = on
    # When On: this will print the anti-stokes average and standard deviation of the peak position for every peak fit
Print_Average_and_STDEV = on
    # When On: if the standard deviation of the peak position if greater than the user input value, a check standard deviation (CHECK STDV FOR) message will be printed along with the file name, sign, and standard deviation of the peak position
Check_STDEV = on
    # When On: this will sort the first column of the original text files
        # If the first row of the first column is less than or equal to zero it will be saved as a sorted file in the Data directory
        # If the first row of the first column is greater than zero it will be saved as a sorted file in the Pressure Calibration directory
Sort_and_Save_Files = on
    # When On: this will analyze the unsorted files in the original directory inputed by the user (Old_Folder)
Analyse_Original_Files = on
# Turn on or off which peak fits (Gaussian, Lorentzian, and Voigt) along with which background fits (Constant, Linear, and Second-Order Polynomial) to analyze
enabled_fits_and_backgrounds = {
    "Gaussian": {
        "Constant": on,
        "Linear": on,
        "Second-Order Polynomial": on
    },
    "Lorentzian": {
        "Constant": on,
        "Linear": on,
        "Second-Order Polynomial": on
    },
    "Voigt": {
        "Constant": on,
        "Linear": on,
        "Second-Order Polynomial": on
    }
}
# Turn on or off the positive or negative peak fits to analyze
enabled_signs = {
    "Positive": on,
    "Negative": on
}
##### ONLY USE POSITIVE NUMBERS #####
# List holding the mask values for the peak fit
    # Lower Mask, Upper Mask
Data_Mask = [[160, center_guess+60], []]
# List holding the mask values for the background fit
    # Lower Mask Before peak, Upper Mask Before peak, Lower Mask After peak, Upper Mask After peak
Background_Mask = [[160, center_guess-20, center_guess+20, 200], []]
# Dictionary holding initial guesses for peaks positive then negative
initial_guess_peak = {
    "Gaussian": [[600, center_guess, 4], [600, -center_guess, 4]],
    "Lorentzian": [[600, center_guess, 5], [600, -center_guess, 5]],
    "Voigt": [[600, center_guess, 4, 5], [600, -center_guess, 4, 5]]
}
# Dictionary holding initial guesses for background
initial_guess_background = {
    "Gaussian": {
        "Constant": [800],
        "Linear": [0.5, 2],
        "Second-Order Polynomial": [0.1, 0.5, 2]
    },
    "Lorentzian": {
        "Constant": [800],
        "Linear": [0.5, 2],
        "Second-Order Polynomial": [0.1, 0.5, 2]
    },
    "Voigt": {
        "Constant": [800],
        "Linear": [0.5, 2],
        "Second-Order Polynomial": [0.1, 0.5, 2]
    }
}

#######################################################################
# End Of User Input Required
#######################################################################

# Dictionary holding positive and negative indicators
Signs = {
    "Positive": 0,
    "Negative": 1
}
# This will make sure the user input for Analyse_Original_Files is followed
if Analyse_Original_Files:
    Analyse_Sorted_Files = off
if not Analyse_Original_Files:
    Analyse_Sorted_Files = on
# Function to fill the negative values of the user inputed initial guesses
def fill_negative_values(data):
    if isinstance(data, dict):
        for key, value in data.items():
            fill_negative_values(value)
    elif isinstance(data, list) and len(data) == 2 and not data[1]:
        data[1] = [-v for v in data[0]]
fill_negative_values(Data_Mask)
fill_negative_values(Background_Mask)
# Create directories for the sorted files and output CSV file and images
dir_list = []
dir_list = os.listdir("%s/%s" % (Path, Old_Folder))
number_of_files = int(len(dir_list))
text_file = "txt"
list_of_text_files = []
list_of_negative_first_number_text_files = []
list_of_positive_first_number_text_files = []
for file in dir_list:
    file_line = file.split(".")
    file_type = file_line[1]
    if file_type == text_file:
        list_of_text_files.append(file)
# For every file in the original directory inputed by the user (Old_Folder) sort the first column of the text file and then print it to a new sorted file
if Sort_and_Save_Files:
    for file in list_of_text_files:
        array = []
        # Load the text file in as an array
        array = loadtxt("%s/%s/%s" % (Path, Old_Folder, file))
        # Create array and sort by the first column
        array = np.array(array)
        array = array[np.lexsort((array[:,0], array[:,0]))][::1]
        # Check if the first row of the first column has a negative number
        if float(array[0,0]) <= 0:
            list_of_negative_first_number_text_files.append(file)
            # Save sorted array a new 'Sorted' file in the Data directory
            if not os.path.exists("%s/%s/Data" % (Path, New_Folder)):
                os.makedirs("%s/%s/Data" % (Path, New_Folder))
            output_file = open("%s/%s/Data/Sorted_%s" % (Path, New_Folder, file), "w")  
            for num in array:
                first_column = num[0]     # x-axis
                second_column = num[1]     # intensity
                output_file.write("%0.16f\t%0.16f\n" % (first_column, second_column))
            output_file.close()
            if not os.path.exists("%s/%s/Data/Images" % (Path, New_Folder)):
                os.makedirs("%s/%s/Data/Images" % (Path, New_Folder))
        else:
            list_of_positive_first_number_text_files.append(file)
            # Save sorted array as a new 'Sorted' file in the Pressure Calibration directory
            if not os.path.exists("%s/%s/Pressure_Calibration" % (Path, New_Folder)):
                os.makedirs("%s/%s/Pressure_Calibration" % (Path, New_Folder))  
            output_file = open("%s/%s/Pressure_Calibration/Sorted_%s" % (Path, New_Folder, file), "w")
            for num in array:
                first_column = num[0]     # x-axis
                second_column = num[1]     # intensity
                output_file.write("%0.16f\t%0.16f\n" % (first_column, second_column))
            output_file.close()        
    output_file.close()
# Create the CSV file and write the header row
file_line = file.split(".")
if not os.path.exists("%s/%s/Data/%s" % (Path, New_Folder, CSV_File_Name)):
    csv_file = open("%s/%s/Data/%s" % (Path, New_Folder, CSV_File_Name), "w", newline='')
    writer = csv.writer(csv_file)
    row = []
    row = ["File Name", 
           "Time Processed", 
           "Average Peak Position Positive", "Standard Deviation of Peak Position Positive", "Average FWHM Positive",
           "Average Peak Position Negative", "Standard Devation of Peak Position Negative", "Average FWHM Negative",
           "Peak Position", "FWHM", "Gaussian", "Constant", 
           "Peak Position", "FWHM", "Gaussian", "Linear", 
           "Peak Position", "FWHM", "Gaussian", "Second-Order Polynomial", 
           "Peak Position", "FWHM", "Lorentzian", "Constant", 
           "Peak Position", "FWHM", "Lorentzian", "Linear", 
           "Peak Position", "FWHM", "Lorentzian", "Second-Order Polynomial", 
           "Peak Position", "FWHM", "Voigt", "Constant", 
           "Peak Position", "FWHM", "Voigt", "Linear", 
           "Peak Position", "FWHM", "Voigt", "Second-Order Polynomial", 
           "Peak Position", "FWHM", "Gaussian", "Constant", 
           "Peak Position", "FWHM", "Gaussian", "Linear", 
           "Peak Position", "FWHM", "Gaussian", "Second-Order Polynomial", 
           "Peak Position", "FWHM", "Lorentzian", "Constant", 
           "Peak Position", "FWHM", "Lorentzian", "Linear", 
           "Peak Position", "FWHM", "Lorentzian", "Second-Order Polynomial", 
           "Peak Position", "FWHM", "Voigt", "Constant", 
           "Peak Position", "FWHM", "Voigt", "Linear", 
           "Peak Position", "FWHM", "Voigt", "Second-Order Polynomial"]
    writer.writerow(row)
    csv_file.close()
# Background Functions
def constant_background(x, c):
    return np.full_like(x, c)
def linear_background(x, m, c):
    return m * x + c
def polynomial_background(x, a, b, c):
    return a * x**2 + b * x + c
background_functions = {
    "Constant": constant_background,
    "Linear": linear_background,
    "Second-Order Polynomial": polynomial_background
}
# Peak Fitting Functions
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
def lorentzian(x, amp, x0, gamma):
    return amp * (gamma**2) / ((x - x0)**2 + gamma**2)
def voigt(x, amp, x0, sigma, gamma):
    z = ((x - x0) + 1j*gamma) / (sigma * np.sqrt(2))
    return amp * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))
peak_fitting_functions = {
    "Gaussian": gaussian,
    "Lorentzian": lorentzian,
    "Voigt": voigt
}
peak_fitting_bounds = {
    # Lower Bound, Upper Bound
    "Gaussian": ([0, -np.inf, 0], [np.inf, np.inf, np.inf]),    # (amp >0), (x0 > -inf), (sigma > 0)
    "Lorentzian": ([0, -np.inf, 0], [np.inf, np.inf, np.inf]),  # (amp >0), (x0 > -inf), (gamma > 0)
    "Voigt": ([0, -np.inf, 0, 0], [np.inf, np.inf, np.inf, np.inf]) # (amp >0), (x0 > -inf), (sigma > 0), (gamma > 0)
}
# Fit Background and Peaks
def fit_peak(x, y, background_name, background_func, peak_name, peak_func, initial_guess_background, initial_guess_peak, sign_name, sign_position, maxfev=500000, max_iterations=500, tolerance=0.01):
    # Fit background
    if sign_name == "Positive":
        mask = (x >= Background_Mask[sign_position][0]) & (x <= Background_Mask[sign_position][1]) | (x >= Background_Mask[sign_position][2]) & (x <= Background_Mask[sign_position][3])
    if sign_name == "Negative":
        mask = (x >= Background_Mask[sign_position][1]) & (x <= Background_Mask[sign_position][0]) | (x >= Background_Mask[sign_position][3]) & (x <= Background_Mask[sign_position][2])
    x_background = x[mask]
    y_background = y[mask]
    popt_background, _ = curve_fit(background_func, x_background, y_background, p0=initial_guess_background, maxfev=maxfev)
    fitted_background = background_func(x, *popt_background)
    y_corrected = y - fitted_background
    bounds = peak_fitting_bounds[peak_name]
    # Iterative fitting for peak
    converged = False
    iterations = 0
    while not converged and iterations < max_iterations:
        iterations += 1
        popt_peak, _ = curve_fit(peak_func, x, y_corrected, p0=initial_guess_peak, bounds=bounds, maxfev=maxfev)
        if np.abs(initial_guess_peak[1] - popt_peak[1]) < tolerance:
            converged = True
        initial_guess_peak = popt_peak
    if not converged and iterations >= max_iterations:
        print(f"{file}: FAIL - {sign_name}: {peak_name} peak Fit: {background_name} Background")
    return popt_peak, y_corrected, fitted_background, iterations, converged
# Plotting a graph for each peak fit and background combination
def plot_fit(x, y, y_corrected, fitted_background, popt_peak, peak_func, peak_name, background_name, file, sign_name):
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.scatter(x, y, color='dimgray', s=5, label='Data')
    plt.plot(x, fitted_background, color='red', label=f"{background_name}")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Relative Wavenumbers (cm-1)')
    plt.ylabel('')
    plt.title(f"Data with {background_name} Background")
    plt.subplot(2, 1, 2)
    plt.scatter(x, y_corrected, color='dimgray', s=5, label='Background-Corrected Data')
    plt.plot(x, peak_func(x, *popt_peak), color='green', label='Fitted Gaussian')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Relative Wavenumbers (cm-1)')
    plt.ylabel('')
    file_line = file.split(".")
    file_line = file_line[0].replace("Sorted_", "")
    plt.title(f"Background-Corrected Data with {peak_name} Peak Fit")
    plt.tight_layout()
    plt.savefig(f"{Path}/{New_Folder}/Data/Images/{file_line}_{peak_name}_Peak_Fit_with_a_{background_name}_Background_{sign_name}", dpi=300, bbox_inches='tight')
# Plotting a comparison graph of all peak fit and background combinations for each text file
def plot_all_fits(file, sign_name, all_fits):
    plt.figure(figsize=(10, 6))
    sign_data = all_fits[file]
    for sign_name_in, lines in sign_data.items():
        if sign_name_in == sign_name:
            if enabled_signs[sign_name]:
                for peak_name, backgrounds in lines.items():
                    for background_name, values in backgrounds.items():
                        if peak_name in enabled_fits_and_backgrounds and background_name in enabled_fits_and_backgrounds[peak_name] and enabled_fits_and_backgrounds[peak_name][background_name]:
                            x = values['x']
                            y = values['y']
                            y_corrected = values['y_corrected']
                            fitted_background = values['fitted_background']
                            popt_peak = values['popt_peak']
                            peak_func = values['peak_func']
                            plt.plot(x, peak_func(x, *popt_peak), label=f'{peak_name} {background_name}')
    plt.scatter(x, y_corrected, color='dimgray', s=5, label='Data')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Relative Wavenumbers (cm-1)')
    plt.ylabel('')
    file_line = file.split(".")
    file_line = file_line[0].replace("Sorted_", "")
    plt.title(f"Peak Fits for {file_line} - {sign_name}")
    plt.tight_layout()
    plt.savefig(f"{Path}/{New_Folder}/Data/Images/{file_line}_Comparison_Plot_{sign_name}", dpi=300, bbox_inches='tight')

#####  
# Main processing function
#####
def process_files(path, old_folder, new_folder, csv_file_name):
    if Analyse_Original_Files:
        data_dir_list = os.listdir(f"{path}/{old_folder}")
    if Analyse_Sorted_Files:
        data_dir_list = os.listdir(f"{path}/{new_folder}/Data")
    all_fits = {}
    for file in data_dir_list:
        if file not in ["Images", csv_file_name]:
            all_fits[file] = {}
            row = []
            statistics_list = []
            peak_info_list = []
            row.append(file)
            if Analyse_Original_Files:
                file_path = f"{path}/{old_folder}/{file}"
            if Analyse_Sorted_Files:
                file_path = f"{path}/{new_folder}/Data/{file}"
            for sign_name, sign_position in Signs.items():
                if enabled_signs[sign_name]:
                    all_fits[file][sign_name] = {}
                    array = np.loadtxt(file_path)
                    x, y = array[:, 0], array[:, 1]
                    if sign_name == "Positive":
                        mask = (x >= Data_Mask[sign_position][0]) & (x <= Data_Mask[sign_position][1])
                    if sign_name == "Negative":
                        mask = (x >= Data_Mask[sign_position][1]) & (x <= Data_Mask[sign_position][0])
                    x, y = x[mask], y[mask]
                    Peak_Position = []
                    FWHM_List = []
                    # Fit peaks with each background function
                    for peak_name, peak_type in peak_fitting_functions.items():
                        all_fits[file][sign_name][peak_name] = {}
                        for background_name, background_function in background_functions.items():
                            if peak_name in enabled_fits_and_backgrounds and background_name in enabled_fits_and_backgrounds[peak_name] and enabled_fits_and_backgrounds[peak_name][background_name]:
                                all_fits[file][sign_name][peak_name][background_name] = {}
                                background_type = background_function.__name__
                                initial_guess_bg = initial_guess_background[peak_name][background_name]
                                initial_guess_crv = initial_guess_peak[peak_name]
                                popt_peak, y_corrected, fitted_background, iterations, converged = fit_peak(x, y, background_name, background_function, peak_name, peak_type, initial_guess_bg, initial_guess_crv, sign_name, sign_position)
                                if peak_name == "Gaussian":
                                    a, x0, sigma = popt_peak
                                    fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma
                                if peak_name == "Lorentzian":
                                    amp, x0, gamma = popt_peak
                                    fwhm = 2 * gamma
                                if peak_name == "Voigt":
                                    amp, x0, sigma, gamma = popt_peak
                                    fwhm = 0.5346 * 2 * gamma + np.sqrt(0.2166 * (2 * gamma)**2 + (2 * sigma * np.sqrt(2 * np.log(2)))**2)
                                peak_center_position = x0
                                full_width_half_max = fwhm
                                FWHM_List.append(fwhm)
                                peak_info_list.append(x0)
                                peak_info_list.append(fwhm)
                                peak_info_list.append(peak_name)
                                peak_info_list.append(background_name)
                                Peak_Position.append(x0)
                                # Store fit results
                                all_fits[file][sign_name][peak_name][background_name] = {
                                    'x': x,
                                    'y': y,
                                    'y_corrected': y_corrected,
                                    'fitted_background': fitted_background,
                                    'popt_peak': popt_peak,
                                    'peak_func': peak_type
                                }
                                # Run the function that will create a graph for each peak fit and background combination
                                if Make_Images:
                                    plot_fit(x, y, y_corrected, fitted_background, popt_peak, peak_type, peak_name, background_name, file, sign_name)
                            else:
                                for _ in range(4):
                                    peak_info_list.append("Skipped")
                else:
                    for _ in range(36):
                        peak_info_list.append("Skipped")
                Average_Peak_Position = np.average(Peak_Position)
                Standard_Deviation_of_Peak_Position = np.std(Peak_Position)
                Average_FWHM = np.average(FWHM_List)
                if enabled_signs[sign_name]:
                    statistics_list.append(Average_Peak_Position)
                    statistics_list.append(Standard_Deviation_of_Peak_Position)
                    statistics_list.append(Average_FWHM)
                    # If the standard deviation of the peak position if greater than the user input value, a check standard deviation (CHECK STDV FOR) message will be printed along with the file name, sign, and standard deviation of the peak position
                    if Standard_Deviation_of_Peak_Position >= Check_Standard_Deviation_of_Peak_Position:
                        if Check_STDEV:
                            print("CHECK STDEV FOR")
                            print(f"\t{file} - {sign_name}")
                            print(f"\tSTDEV = {Standard_Deviation_of_Peak_Position:0.5f}")
                            print()
                    # This will print the anti-stokes average and standard deviation of the peak position for every peak fit
                    if Print_Average_and_STDEV:
                        print(f"{file} - {sign_name}")
                        print(f"Anti-Stokes Average = {Average_Peak_Position:0.5f} - STDEV = {Standard_Deviation_of_Peak_Position:0.5f}")
                        print()
                    # Run the function that will create a comparison graph of all peak fit and background combinations for each text file
                    if Make_Comparison_Graph:
                        plot_all_fits(file, sign_name, all_fits)
                else:
                    for _ in range(3):
                        statistics_list.append("Skipped")
            # Write the cure fits for this file to the CSV file
            time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            row.append(time)
            row = row + statistics_list + peak_info_list
            with open(f"{Path}/{New_Folder}/Data/{CSV_File_Name}", "a", newline='') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(row)
#end

# Run the main processing function
process_files(Path, Old_Folder, New_Folder, CSV_File_Name)


# Copyright (C) 2025 Pease et al.
# Last Updated: March 25th, 2025


# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
