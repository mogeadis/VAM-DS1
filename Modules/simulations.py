# ====================================================================================================
# # ==================================================================================================
# # Filename: simulations.py
# #
# # Description: This module contains functions and classes used for simulating circuits
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# Import Modules (Built-In)
import os
from abc import ABC,abstractmethod
from itertools import product
from collections.abc import Iterable

# Import Modules (User-Defined)
import utilities as utils

# Import Modules (Third-Party)
import numpy as np
import pandas as pd
import scipy.fft as fft
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,ScalarFormatter

# ====================================================================================================
# Function: frequencySpectrum()
# Description: This function computes the frequency spectrum of a signal
# ====================================================================================================
def frequencySpectrum(signal,fs):

    # Compute Complex Spectrum
    spectrum = fft.fft(signal)
    frequencies = fft.fftfreq(len(signal),1/fs)

    # Keep Real Frequencies
    freqs = frequencies[:len(frequencies)//2]
    spectrum = spectrum[:len(frequencies)//2]

    # Calculate Magnitude & Phase
    magnitude = 20*np.log10(np.abs(spectrum))
    phase = unwrapPhase(np.degrees(np.angle(spectrum)))

    # Return
    return magnitude,phase,freqs

# ====================================================================================================
# Function: unwrapPhase()
# Description: This function unwraps the phase response of a signal
# ====================================================================================================
def unwrapPhase(phase):

    # Compute Unwrapped Phase
    phase = np.unwrap(phase,period = 360)
    if phase[1] >= 180:
        phase -= 360
    elif phase[1] <= -180:
        phase += 360

    # Return
    return phase

# ====================================================================================================
# Class: Trace
# Description: This class represents a 2D plot trace with data and style information
# ====================================================================================================
class Trace():

    # Constructor Method
    def __init__(self,x,y,linestyle = '-',linewidth = 2,color = 'k',label = '',tag = ''):

        # Assign Arguments
        self.x = x
        self.y = y
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.color = color
        self.label = label
        self.tag = tag

# ====================================================================================================
# Class: BaseSource
# Description: This class serves as a base for signal source classes
# ====================================================================================================
class BaseSource(ABC):

    # Constructor Method
    @abstractmethod
    def __init__(self,sample_rate = 44100,oversampling = 1):

        # Assign Arguments
        self.sample_rate = sample_rate
        self.oversampling = oversampling

        # Define Member Variables
        self.fs = self.sample_rate*self.oversampling

    # Public Methods
    @abstractmethod
    def generate(self,simulation_time):
        pass

    # Private Methods
    def createLabel(self):
        return ''

# ====================================================================================================
# Class: SineSource
# Description: This class generates a sine wave source
# ====================================================================================================
class SineSource(BaseSource):

    # Constructor Method
    def __init__(self,amplitude = 1,frequency = 1000,phase = 0,dc_offset = 0,sample_rate = 44100,oversampling = 1):
        
        # Initialize Base Class
        super().__init__(sample_rate,oversampling)

        # Assign Arguments
        self.amplitude = amplitude
        self.frequency = frequency
        self.phase = phase
        self.dc_offset = dc_offset

        # Define Member Variables
        self.label = self.createLabel()

    # Public Methods
    def generate(self,simulation_time):

        # Generate Sine Wave
        t = np.arange(0,simulation_time,1/self.fs)
        signal = self.dc_offset + self.amplitude*np.sin(2*np.pi*self.frequency*t + self.phase)

        # Return
        return signal,t
    
    # Private Methods
    def createLabel(self):

        # Create Label String
        label = ''
        if self.dc_offset != 0:
            label += f'{self.dc_offset} + '
        if self.amplitude != 1:
            label += f'{self.amplitude}'
        label += 'sin('
        if self.frequency != 0:
            label += f'2π{self.frequency}t'
            if self.phase != 0:
                label += f' + {self.phase}'
        else:
            label += f'{self.phase}'
        label += ')'

        # Return
        return label
        
# ====================================================================================================
# Class: ImpulseSource
# Description: This class generates a unit impulse source
# ====================================================================================================
class ImpulseSource(BaseSource):

    # Constructor Method
    def __init__(self,sample_rate,oversampling = 1,amplitude = 1):

        # Initialize Base Class
        super().__init__(sample_rate,oversampling)

        # Assign Arguments
        self.amplitude = amplitude

        # Define Member Variables
        self.label = self.createLabel()

    # Public Methods
    def generate(self,simulation_time):

        # Generate Unit Impulse
        t = np.arange(0,simulation_time,1/self.fs)
        signal = np.array([self.amplitude] + (len(t) - 1)*[0])

        # Return
        return signal,t
    
    # Private Methods
    def createLabel(self):

        # Create Label String
        label = ''
        if self.amplitude != 1:
            label += f'{self.amplitude}'
        label += 'δ(t)'

        # Return
        return label

# ====================================================================================================
# Class: BaseAnalysis
# Description: This class serves as a base for circuit analysis classes
# ====================================================================================================
class BaseAnalysis(ABC):

    # Static Variables
    @property
    def analysis_type():
        return ''

    # Constructor Method
    def __init__(self,model,source,parameters,sample_rate,oversampling = 1,simulation_time = 1,spice_simulations = False,linearize = False,):
        
        # Assign Arguments
        self.model = model
        self.source = source
        self.parameters = parameters
        self.sample_rate = sample_rate
        self.oversampling = oversampling
        self.simulation_time = simulation_time
        self.spice_simulations = spice_simulations
        self.linearize = linearize
        
        # Define Member Variables
        self.fs = self.sample_rate*self.oversampling
        self.spice_info = self.fetchSimulationInfo()
        self.steps = self.generateSteps()
        self.traces = []

    # Public Methods
    @abstractmethod
    def run(self):
        pass

    def plot(self):
        pass

    # Private Methods
    def createLabel(self,step):

        # Create Label String
        label = ''
        for index,(key,value) in enumerate(step.items()):
            if isinstance(self.parameters[key],Iterable):
                label += f'{key.capitalize()} = {int(100*value)}% | '

        # Return
        return label
    
    def createTitle(self):

        # Create Title String
        title = f'{self.model.name} | {self.__class__.analysis_type} Analysis\n'
        title += f'Sample Rate: {self.sample_rate/1000}kHz | Oversampling: $\\times${self.oversampling}'
        for key,value in self.parameters.items():
            if not isinstance(value,Iterable):
                title += f' | {key.capitalize()}: {int(100*value)}%'

        # Return
        return title
    
    def generateSteps(self):

        # Initialize Variables
        fixed_params = {}
        step_param_keys = []
        step_param_values = []

        # Analyze Parameters
        for key,value in self.parameters.items():
            if isinstance(value,Iterable):
                step_param_keys.append(key)
                step_param_values.append(value)
            else:
                fixed_params[key] = value

        # Generate Parameter Steps
        combos = list(product(*step_param_values))
        steps = [{**dict(zip(step_param_keys,combo)),**fixed_params} for combo in combos]

        # Return
        return steps
    
    def fetchSimulationInfo(self):

        # Check SPICE Flag
        if self.spice_simulations:

            # Load SPICE Simulations Info
            spice_info = pd.read_csv(utils.INFO_CSV)

            # Execute Query
            prefix = self.model.name.replace(' ','').replace('-','') + f'_{self.__class__.analysis_type[0:2].upper()}_'
            string = f'ID.str.startswith(\'{prefix}\')'
            if isinstance(self.source,SineSource):
                string += f' and AMPLITUDE == {self.source.amplitude} and FREQUENCY == {self.source.frequency}'
            spice_info = spice_info.query(string)

            # Check Integrity
            if spice_info.empty:
                print('No available SPICE simulations')
                self.spice_simulations = False
        else:
            spice_info = pd.DataFrame()

        return spice_info

# ====================================================================================================
# Class: TransientAnalysis
# Description: This class simulates a Transient Analysis for a given circuit
# ====================================================================================================
class TransientAnalysis(BaseAnalysis):

    # Static Variables
    analysis_type = 'Transient'

    # Constructor Class
    def __init__(self,model,source,parameters,sample_rate,oversampling = 1,simulation_time = 1,spice_simulations = True,linearize = False):

        # Initialize Base Class
        super().__init__(**{key : value for key,value in locals().items() if key not in ["self","__class__"]})

    # Public Methods
    def run(self):

        # Generate Input Signal
        Vin,t = self.source.generate(self.simulation_time)
        N = len(Vin)
        self.traces.append(Trace(x = t,y = Vin,linestyle = '-',linewidth = 3.5,color = 'k',label = f'$V_{{in}}$ = {self.source.label}'))

        # Iterate Over Steps
        for index,step in enumerate(self.steps):

            # Run Simulation
            print(f'Step {index + 1}/{len(self.steps)}: Running Simulation ...',end = ' ',flush = True)
            model = self.model(fs = self.fs,linearize = self.linearize,**step)
            Vout = np.zeros(N)
            for n in range(N):
                Vout[n] = model.processSample(Vin[n])

            # Store Simulation Data
            step_label = self.createLabel(step)
            step_color = plt.get_cmap('tab10')(index % 10)
            self.traces.append(Trace(x = t,y = Vout,linestyle = '-',linewidth = 2.5,color = step_color,label = step_label + f'{model.label}'))

            # Check SPICE Flag
            if self.spice_simulations:

                # Execute Query
                query_string = ''
                for index,(key,value) in enumerate(step.items()):
                    query_string += f'{key.upper()} == {value}'
                    if index < len(step) - 1:
                        query_string += ' and '
                if query_string:
                    step_info = self.spice_info.query(query_string)
                else:
                    step_info = self.spice_info

                # Check Integrity
                if not step_info.empty:

                    # Load SPICE Data
                    simulation_id = step_info.iloc[0]['ID']
                    spice_data = pd.read_csv(os.path.join(utils.SIMULATIONS_CSV_PATH,simulation_id + '.csv'))
                    self.traces.append(Trace(x = spice_data['TIME'],y = spice_data['AMPLITUDE'],linestyle = '--',linewidth = 2,color = utils.adjustLightness(step_color),label = step_label + 'SPICE'))
                else:
                    print(f'No available SPICE simulations ...',end = ' ',flush = True)
            print('Analysis Completed!')

    def plot(self,scaling = 1,x_unit = 's',x_min = None,x_max = None,y_unit = 'V',y_min = None,y_max = None,title = None):

        # Configure Axis Units
        x_factor = utils.unitPrefixFactor(x_unit)
        y_factor = utils.unitPrefixFactor(y_unit)

        # Plot Data
        fig,ax = plt.subplots(figsize = (6.4*scaling,4.8*scaling))
        for trace in self.traces:
            ax.plot(trace.x*x_factor,trace.y*y_factor,color = trace.color,linestyle = trace.linestyle,linewidth = trace.linewidth,label = trace.label)

        # Adjust x-Axis Limits
        if x_min is None:
            x_min = 0
        if x_max is None:
            x_max = self.simulation_time*x_factor
        ax.set_xlim(x_min,x_max)

        # Adjust y-Axis Limits
        if y_min is not None:
            ax.set_ylim(bottom = y_min)
        if y_max is not None:
            ax.set_ylim(top = y_max)

        # Set Axis Labels
        ax.set_xlabel(f'Time [{x_unit}]')
        ax.set_ylabel(f'Amplitude [{y_unit}]')

        # Enable Grid & Legend
        ax.grid(visible = True)
        ax.legend(loc = 'lower right',framealpha = 0.9)

        # Set Title
        if title is None:
            title = self.createTitle()
        ax.set_title(title)
        fig.canvas.manager.set_window_title(title.replace('\n',' | ').replace(r'$\times$','x'))

# ====================================================================================================
# Class: ACAnalysis
# Description: This class simulates an AC Analysis for a given circuit
# ====================================================================================================
class ACAnalysis(BaseAnalysis):

    # Static Variables
    analysis_type = 'AC'

    # Constructor Method
    def __init__(self,model,source,parameters,sample_rate,oversampling = 1,simulation_time = 1,spice_simulations = True,linearize = False):

        # Initialize Base Class
        super().__init__(**{key : value for key,value in locals().items() if key not in ["self","__class__"]})

    # Public Methods
    def run(self):

        # Generate Input Signal
        Vin = self.source.generate(self.simulation_time)[0]
        N = len(Vin)

        # Iterate Over Steps
        for index,step in enumerate(self.steps):

            # Run Simulation
            print(f'Step {index + 1}/{len(self.steps)}: Running Simulation ...',end = ' ',flush = True)
            model = self.model(fs = self.fs,linearize = self.linearize,**step)
            Vout = np.zeros(N)
            for n in range(N):
                Vout[n] = model.processSample(Vin[n])

            # Compute Frequency Response
            magnitude,phase,freqs = frequencySpectrum(Vout,self.fs)

            # Store Simulation Data
            step_label = self.createLabel(step)
            step_color = plt.get_cmap('tab10')(index % 10)
            self.traces.append(Trace(x = freqs,y = magnitude,linestyle = '-',linewidth = 2.5,color = step_color,label = step_label + f'{model.label}',tag = 'magnitude'))
            self.traces.append(Trace(x = freqs,y = phase,linestyle = '-',linewidth = 2.5,color = step_color,label = step_label + f'{model.label}',tag = 'phase'))

            # Check SPICE Flag
            if self.spice_simulations:

                # Execute Query
                query_string = ''
                for index,(key,value) in enumerate(step.items()):
                    query_string += f'{key.upper()} == {value}'
                    if index < len(step) - 1:
                        query_string += ' and '
                if query_string:
                    step_info = self.spice_info.query(query_string)
                else:
                    step_info = self.spice_info

                # Check Integrity
                if not step_info.empty:

                    # Load SPICE Data
                    simulation_id = step_info.iloc[0]['ID']
                    spice_data = pd.read_csv(os.path.join(utils.SIMULATIONS_CSV_PATH,simulation_id + '.csv'))
                    self.traces.append(Trace(x = spice_data['FREQUENCY'],y = spice_data['MAGNITUDE'],linestyle = '--',linewidth = 2,color = utils.adjustLightness(step_color),label = step_label + 'SPICE',tag = 'magnitude'))
                    self.traces.append(Trace(x = spice_data['FREQUENCY'],y = unwrapPhase(spice_data['PHASE']),linestyle = '--',linewidth = 2,color = utils.adjustLightness(step_color),label = step_label + 'SPICE',tag = 'phase'))
                else:
                    print(f'No available SPICE simulations ...',end = ' ',flush = True)
            print('Analysis Completed!')

    def plot(self,scaling = 1,x_unit = 'Hz',x_min = None,x_max = None,y_unit = 'dB',y_min = None,y_max = None,title = None):

        # Configure x-Axis Unit
        if x_unit == 'Hz':
            x_factor = 1
        elif x_unit == 'rad/s':
            x_factor = 2*np.pi/(self.sample_rate)
        else:
            raise ValueError('Invalid unit (x axis)')

        # Configure y-Axis Unit
        if y_unit == 'dB':
            plot_type = 'magnitude'
            y_factor = 1
        elif y_unit == 'degrees':
            plot_type = 'phase'
            y_factor = 1
            y_unit = '\N{DEGREE SIGN}'
        elif y_unit == 'rad':
            plot_type = 'phase'
            y_factor = np.pi/180
        else:
            raise ValueError('Invalid unit (y axis)')

        # Plot Data
        fig,ax = plt.subplots(figsize = (6.4*scaling,4.8*scaling))
        for trace in self.traces:
            if trace.tag == plot_type:
                ax.plot(trace.x*x_factor,trace.y*y_factor,color = trace.color,linestyle = trace.linestyle,linewidth = trace.linewidth,label = trace.label)

        # Adjust x-Axis Limits
        ax.set_xscale('log')
        if x_unit == 'Hz':
            if x_min is None:
                x_min = 10
            if x_max is None:
                x_max = self.sample_rate/2
            ax.xaxis.set_major_formatter(FuncFormatter(lambda f,_: f'{int(f)}' if f < 1000 else f'{int(f/1000)}k'))
        elif x_unit == 'rad/s':
            if x_min is None:
                x_min = 10*x_factor
            if x_max is None:
                x_max = np.pi
            ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xlim(x_min,x_max)

        # Adjust y-Axis Limits
        if y_min is not None:
            ax.set_ylim(bottom = y_min)
        if y_max is not None:
            ax.set_ylim(top = y_max)

        # Set Axis Labels
        ax.set_xlabel(f'Frequency [{x_unit}]')
        ax.set_ylabel(f'{plot_type.capitalize()} [{y_unit}]')

        # Enable Grid & Legend
        ax.grid(visible = True,which = 'major')
        ax.grid(visible = True,which = 'minor',alpha = 0.5)
        ax.legend(loc = 'best',framealpha = 0.9)

        # Set Title
        if title is None:
            title = self.createTitle().replace('\n',f' | {plot_type.capitalize()} Response\n')
        ax.set_title(title)
        fig.canvas.manager.set_window_title(title.replace('\n',' | ').replace(r'$\times$','x'))