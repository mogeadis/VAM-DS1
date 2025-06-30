# ====================================================================================================
# # ==================================================================================================
# # Filename: data_formatting.py
# #
# # Description: This script formats the data exported by LTspice into a .csv file
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# ====================================================================================================
# Execution Initialization
# ====================================================================================================

# Import Modules (Built-In)
import os
import sys
import time

# Import Modules (User-Defined)
sys.path.append(os.path.abspath('Modules'))
import utilities as utils

# Import Modules (Third-Party)
None

# Start Timer
start_time = time.perf_counter()
utils.printHeader('data_formatting.py')

# ====================================================================================================
# Data Formatting
# ====================================================================================================

# Access Simulation Files
data_directory = os.listdir(utils.SIMULATIONS_TXT_PATH)
print('Data Directory:')
print(utils.SIMULATIONS_TXT_PATH,'\n')
for file_name in data_directory:

    # Check File Type
    file_stem,file_extension = os.path.splitext(file_name)
    if file_extension == '.txt':

        # Create Data Header
        data = ''
        if '_TR_' in file_stem:
            data += 'TIME,AMPLITUDE\n'
        elif '_AC_' in file_stem:
            data += 'FREQUENCY,MAGNITUDE,PHASE\n'
        else:
            continue

        # Read Data from .txt
        print(f'Formatting \'{file_name}\' ...',end = ' ',flush = True)
        read_path = os.path.join(utils.SIMULATIONS_TXT_PATH,file_name)
        with open(read_path,'r') as read_file:

            # Parse Data
            next(read_file)
            for text_line in read_file:
                text_line = text_line.replace('\t',',')
                text_line = ''.join(([char for char in text_line if char.isdigit() or char in ['e','+','-',',','.']]))
                data += text_line + '\n'

            # Write Data to .csv
            write_path = os.path.join(utils.SIMULATIONS_CSV_PATH,file_stem + '.csv')
            with open(write_path,'w') as write_file:
                write_file.write(data)
                print('Process Completed!')

# ====================================================================================================
# Execution Conclusion
# ====================================================================================================

# End Timer
runtime = utils.calculateElapsedTime(start_time,unit = 'ms')
utils.printFooter(runtime)