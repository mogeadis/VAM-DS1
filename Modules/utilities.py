# ====================================================================================================
# # ==================================================================================================
# # Filename: utilities.py
# #
# # Description: This module contains variables and functions used throughout the project
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# Import Modules (Built-In)
import os
import time
import colorsys
from datetime import datetime

# Import Modules (User-Defined)
None

# Import Modules (Third-Party)
import matplotlib.pyplot as plt

# ====================================================================================================
# Configuration Variables
# ====================================================================================================

# Filesystem Paths
FILES_PATH = os.path.join(os.getcwd(),'Files')
SIMULATIONS_PATH = os.path.join(FILES_PATH,'Simulations')
SIMULATIONS_TXT_PATH = os.path.join(SIMULATIONS_PATH,'TXT')
SIMULATIONS_CSV_PATH = os.path.join(SIMULATIONS_PATH,'CSV')
INFO_CSV = os.path.join(SIMULATIONS_PATH,'info.csv')

# Newton-Raphson Parameters
ITERATIONS = 10
SUBITERATIONS = 10
THRESHOLD = 1e-6

# Circuit Variables
POT_MARGIN = 1e-3
THERMAL_VOLTAGE = 25.85e-3

# ====================================================================================================
# Function: limit()
# Description: This function limits a value between a lower and an upper bound
# ====================================================================================================
def limit(value,lower,upper,margin = 0):

    # Limit Value
    if value <= lower:
        value = lower + margin
    elif value >= upper:
        value = upper - margin

    # Return
    return value

# ====================================================================================================
# Function: adjustLightness()
# Description: This function adjusts the lightness of a color
# ====================================================================================================
def adjustLightness(color,factor = 0.7):

    # Convert RGB to HLS
    h,l,s = colorsys.rgb_to_hls(*color[:3])

    # Adjust Lightness
    l = max(0,min(1,l*factor))

    # Convert HLS to RGB
    color = colorsys.hls_to_rgb(h,l,s) + (color[3],)

    # Return
    return color

# ====================================================================================================
# Function: showPlots()
# Description: This function displays all open figures
# ====================================================================================================
def showPlots(close_figs = True):

    # Draw Figures
    plt.draw() 

    # Pause Execution
    plt.pause(1e-3)
    input('Press enter to continue...\n') 

    # Close Figures
    if close_figs == True:
        plt.close('all')

# ====================================================================================================
# Function: unitPrefixFactor()
# Description: This function returns the prefix and multiplication factor for a given measurement unit
# ====================================================================================================
def unitPrefixFactor(unit):

    # Extract Unit Prefix
    L = len(unit)
    if L > 1:
        prefix = unit[0]
    elif L == 1:
        prefix = ''
    else:
        raise ValueError("Invalid argument")

    # Determine Multiplication Factor
    prefix_factor_dict = {'M' : 1e-6,
                          'k' : 1e-3,
                          'm' : 1e3,
                          'u' : 1e6,
                          'n' : 1e9,
                          'p' : 1e12,
                          'f' : 1e15}
    factor = prefix_factor_dict.get(prefix,1)

    # Return
    return factor

# ====================================================================================================
# Function: calculateElapsedTime()
# Description: This function calculates the elapsed time since a given starting point
# ====================================================================================================
def calculateElapsedTime(start_time,unit = 's',decimals = 3):

    # Adjust Time Unit
    factor = unitPrefixFactor(unit)
    
    # Calculate Elapsed Time
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    elapsed_time = f"{(factor*(elapsed_time)):.{decimals}f} {unit}"

    # Return
    return elapsed_time

# ====================================================================================================
# Function: printTextWithBorders()
# Description: This function prints an aligned text between two border lines
# ====================================================================================================
def printTextWithBorders(text,alignment = 'center',border_char = "=",border_length = 100,line_break = True):

    # Determine Indentation
    border_length = max(border_length,len(text))
    if alignment == 'left':
        indent_size = 0
    elif alignment == 'center':
        indent_size = (border_length - len(text))//2
    elif alignment == 'right':
        indent_size = (border_length - len(text))
    else:
        raise ValueError("Invalid alignment type")
    
    # Print Information
    print("") if line_break else None
    print(border_length*border_char)
    print(indent_size*" " + f"{text}")
    print(border_length*border_char)
    print("") if line_break else None

# ====================================================================================================
# Function: printHeader()
# Description: This function prints a header text, typically at the start of a script execution
# ====================================================================================================
def printHeader(file_name,timestamp = True,**kwargs):

    # Format Text
    text = "Executing File: " + file_name
    text = text + " | " + datetime.now().strftime("%d/%m/%Y %H:%M:%S") if timestamp else text
    printTextWithBorders(text,**kwargs)

# ====================================================================================================
# Function: printFooter()
# Description: This function prints a footer text, typically at the end of a script execution
# ====================================================================================================
def printFooter(runtime,timestamp = True,**kwargs):

    # Format Text
    text = "Execution Runtime: " + runtime
    text = text + " | " + datetime.now().strftime("%d/%m/%Y %H:%M:%S") if timestamp else text
    printTextWithBorders(text,**kwargs)