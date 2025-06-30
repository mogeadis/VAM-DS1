# ====================================================================================================
# # ==================================================================================================
# # Filename: circuit_analysis.py
# #
# # Description: This script simulates a circuit by analyzing its time and frequency response
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
from bossds1 import *
import utilities as utils
import simulations as sim

# Import Modules (Third-Party)
None

# Start Timer
start_time = time.perf_counter()
utils.printHeader('circuit_analysis.py')

# ====================================================================================================
# Circuit Analysis
# ====================================================================================================

# Select Circuit Model
circuit_model = [BossDS1,TransistorAmplifier,OperationalAmplifier,DiodeClipper,PassiveFilter][0]

# Select Analysis Type
analysis_type = ['Transient','AC'][0]

# Set Circuit Parameters
distortion = [0,0.25,0.5,0.75,1][:]
tone = [0,0.25,0.5,0.75,1][2]
level = 1
parameters = {'distortion' : distortion,
                    'tone' : tone,
                   'level' : level}
parameters = {key : parameters[key] for key in circuit_model.parameters}

# Set Sample Rate
sample_rate = 44100
oversampling = 1

# Generate Sine Source
amplitude = [0.2,1][0]
frequency = [100,200,1000,2000][2]
sine_source = sim.SineSource(amplitude = amplitude,frequency = frequency,sample_rate = sample_rate,oversampling = oversampling)

# Generate Impulse Source
impulse_source = sim.ImpulseSource(sample_rate = sample_rate,oversampling = oversampling)

# Run Analysis
spice_simulations = True
if analysis_type == 'Transient':
    print(f'{circuit_model.name} | Transient Analysis:')
    simulation_time = 5/frequency
    linearize = False
    analysis = sim.TransientAnalysis(source = sine_source,model = circuit_model,parameters = parameters,sample_rate = sample_rate,oversampling = oversampling,simulation_time = simulation_time,spice_simulations = spice_simulations,linearize = linearize)
    analysis.run()
    analysis.plot(scaling = 1.3,x_unit = 'ms')
elif analysis_type == 'AC':
    print(f'{circuit_model.name} | AC Analysis:')
    simulation_time = 1
    linearize = True
    analysis = sim.ACAnalysis(source = impulse_source,model = circuit_model,parameters = parameters,sample_rate = sample_rate,oversampling = oversampling,simulation_time = simulation_time,spice_simulations = spice_simulations,linearize = linearize)
    analysis.run()
    analysis.plot(scaling = 1.3,y_unit = 'dB',x_min = 1,x_max = 20000)
    analysis.plot(scaling = 1.3,y_unit = 'degrees',x_min = 1,x_max = 20000)

# ====================================================================================================
# Execution Conclusion
# ====================================================================================================

# End Timer
runtime = utils.calculateElapsedTime(start_time,unit = 's')
utils.printFooter(runtime)

# Display Figures
utils.showPlots()