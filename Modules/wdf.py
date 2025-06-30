# ====================================================================================================
# # ==================================================================================================
# # Filename: wdf.py
# #
# # Description: This module contains classes representaing basic Wave Digital Filter elements
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# Import Modules (Built-In)
from abc import ABC,abstractmethod

# Import Modules (User-Defined)
None

# Import Modules (Third-Party)
import numpy as np

# ====================================================================================================
# Class: BasePort
# Description: This class serves as a base for Wave Digital Filter port elements
# ====================================================================================================
class BasePort(ABC):

    # Constructor Method
    @abstractmethod
    def __init__(self,Rp = 1):

        # Assign Arguments
        self.Rp = Rp

        # Define Member Variables
        self.a = 0
        self.b = 0
        
    # Public Methods
    @abstractmethod
    def getReflectedWave(self):
        pass

    def setIncidentWave(self,a):
        self.a = a

    def setPortResistance(self,Rp):
        self.Rp = Rp

    def getPortResistance(self):
        return self.Rp

    def waveToVoltage(self):

        # Calculate Voltage
        voltage = (self.a + self.b)/2

        # Return
        return voltage
    
    def waveToCurrent(self):

        # Calculate Current
        current = (self.a - self.b)/(2*self.Rp)

        # Return
        return current

# ====================================================================================================
# Class: Resistor
# Description: This class represents a resistor
# ====================================================================================================
class Resistor(BasePort):

    # Constructor Method
    def __init__(self,R):

        # Initialize Base Class
        super().__init__(R)

    # Public Methods
    def getReflectedWave(self):
        self.b = 0
        return self.b

# ====================================================================================================
# Class: Capacitor
# Description: This class represents a capacitor
# ====================================================================================================
class Capacitor(BasePort):

    # Constructor Method
    def __init__(self,C,fs = 44100):

        # Initialize Base Class
        super().__init__(1/(2*fs*C))

    # Public Methods
    def getReflectedWave(self):
        self.b = self.a
        return self.b
    
# ====================================================================================================
# Class: ResistiveVoltageSource
# Description: This class represents a resistive voltage source
# ====================================================================================================
class ResistiveVoltageSource(BasePort):

    # Constructor Method
    def __init__(self,Rs):

        # Initialize Superclass
        super().__init__(Rs)

        # Define Member Variables
        self.Vs = 0

    # Public Methods
    def getReflectedWave(self):
        self.b = self.Vs
        return self.b

    def setVoltage(self,Vs):
        self.Vs = Vs
    
# ====================================================================================================
# Class: ResistiveCurrentSource
# Description: This class represents a resistive current source
# ====================================================================================================
class ResistiveCurrentSource(BasePort):

    # Contructor Method
    def __init__(self,Rs):

        # Initialize Base Class
        super().__init__(Rs)

        # Define Member Variables
        self.Is = 0

    # Public Methods
    def getReflectedWave(self):
        self.b = self.Is*self.Rp
        return self.b
    
    def setCurrent(self,Is):
        self.Is = Is

# ====================================================================================================
# Class: BaseAdaptorRoot
# Description: This class serves as a base for root adaptor elements
# ====================================================================================================
class BaseAdaptorRoot(ABC):

    # Constructor Method
    @abstractmethod
    def __init__(self,*ports):

        # Assign Arguments
        self.ports = ports

        # Define Member Variables
        self.N = len(self.ports)
        self.a = np.zeros(self.N)
        self.b = np.zeros(self.N)
        self.scattering_matrix = np.zeros((self.N,self.N))

    # Public Methods
    @abstractmethod
    def computeScatteringMatrix(self):
        pass

    def collectIncidentWaves(self):
        for index,port in enumerate(self.ports):
            self.a[index] = port.getReflectedWave()

    def propagateReflectedWaves(self):
        self.b = self.scattering_matrix @ self.a
        for index,port in enumerate(self.ports):
            port.setIncidentWave(self.b[index])

# ====================================================================================================
# Class: SeriesAdaptorRoot
# Description: This class represents a root series adaptor
# ====================================================================================================
class SeriesAdaptorRoot(BaseAdaptorRoot):

    # Constructor Method
    def __init__(self,*ports):

        # Initialize Base Class
        super().__init__(*ports)
        self.computeScatteringMatrix()

    # Public Methods
    def computeScatteringMatrix(self):

        # Precompute Sum
        sum = 0
        for port in self.ports:
            sum += port.getPortResistance()

        # Fill Scattering Matrix
        for row in range(self.N):
            coeff = -2*self.ports[row].getPortResistance()/sum
            for col in range(self.N):
                self.scattering_matrix[row, col] = coeff + 1 if row == col else coeff

# ====================================================================================================
# Class: ParallelAdaptorRoot
# Description: This class represents a root parallel adaptor
# ====================================================================================================
class ParallelAdaptorRoot(BaseAdaptorRoot):

    # Constructor Method
    def __init__(self,*ports):

        # Initialize Base Class
        super().__init__(*ports)
        self.computeScatteringMatrix()

    # Public Methods
    def computeScatteringMatrix(self):

        # Precompute Sum
        sum = 0
        for port in self.ports:
            sum += 1/port.getPortResistance()

        # Fill Scattering Matrix
        for col in range(self.N):
            coeff = 2*(1/self.ports[col].getPortResistance())/sum
            for row in range(self.N):
                self.scattering_matrix[row, col] = coeff - 1 if row == col else coeff