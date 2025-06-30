# ====================================================================================================
# # ==================================================================================================
# # Filename: bossds1.py
# #
# # Description: This module contains classes representing models of the Boss DS-1 circuit stages
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# Import Modules (Built-In)
from abc import ABC,abstractmethod

# Import Modules (User-Defined)
import wdf
import nonlinear as nl
import utilities as utils

# Import Modules (Third-Party)
import numpy as np

# ====================================================================================================
# Class: BaseModel
# Description: This class serves as a base for circuit model classes
# ====================================================================================================
class BaseModel(ABC):

    # Static Variables
    @property
    def name(self):
        return self.__class__.__name__
    
    @property
    def label(self):
        return self.__class__.__name__
    
    @property
    def parameters(self):
        return []

    # Constructor Method
    @abstractmethod
    def __init__(self,fs = 44100,**kwargs):
        
        # Assign Arguments
        self.fs = fs

        # Define Member Variables
        self.Ts = 1/self.fs
        self.c = 2/self.Ts

    # Public Methods
    @abstractmethod
    def processSample(self,Vin):
        pass

# ====================================================================================================
# Class: TransistorAmplifier
# Description: This class represents the Boss DS-1 Transistor Amplifier stage
# ====================================================================================================
class TransistorAmplifier(BaseModel):

    # Static Variables
    name = 'Transistor Amplifier'
    label = 'SSM'
    parameters = []

    # Constructor Method
    def __init__(self,fs = 44100,**kwargs):

        # Initialize Base Class
        super().__init__(fs)

        # Dimension Variables
        self.XDIM = 3
        self.VDIM = 2

        # Sources
        self.Vcc = 9

        # Resistors
        self.R6 = 100e3
        self.R7 = 470e3
        self.R8 = 10e3
        self.R9 = 22
        self.R10 = 100e3

        # Conductances
        self.G6 = 1/self.R6
        self.G7 = 1/self.R7
        self.G8 = 1/self.R8
        self.G9 = 1/self.R9
        self.G10 = 1/self.R10

        # Capacitors
        self.C3 = 47e-9
        self.C4 = 250e-12
        self.C5 = 68e-9

        # Transistor
        self.Bf = 200
        self.Br = 0.1
        self.Is = 6.734e-15
        self.n = 1
        self.Vt = utils.THERMAL_VOLTAGE
        self.transistor = nl.Transistor(Bf = self.Bf,Br = self.Br,Is = self.Is,n = self.n,Vt = self.Vt)
        self.linearize = kwargs.get('linearize',False)

        # Newton-Raphson Parameters
        self.iterations = utils.ITERATIONS
        self.subiterations = utils.SUBITERATIONS
        self.threshold = utils.THRESHOLD

        # System Vectors
        self.x = np.zeros((self.XDIM,1))
        self.u = np.array([[0.0,self.Vcc]]).T
        self.v = np.zeros((self.VDIM,1))
        self.fv = np.zeros((self.VDIM,1))
        self.p = np.zeros((self.VDIM,1))

        # System Matrices
        self.A = np.array([[-(self.G6 + self.G8 + self.G10)/self.C3,-(self.G8 + self.G10)/self.C3,-self.G10/self.C3],
                           [-(self.G8 + self.G10)/self.C4,-(self.G7 + self.G8 + self.G10)/self.C4,-self.G10/self.C4],
                           [-self.G10/self.C5,-self.G10/self.C5,-self.G10/self.C5]])
        self.B = np.array([[(self.G6 + self.G8 + self.G10)/self.C3,-self.G8/self.C3],
                           [(self.G8 + self.G10)/self.C4,-self.G8/self.C4],
                           [self.G10/self.C5,0]])
        self.C = np.array([[1/self.C3,1/self.C3],
                           [0,1/self.C4],
                           [0,0]])
        self.D = np.array([[-1,0,0],
                           [0,1,0]])
        self.E = np.array([[1,0],
                           [0,0]])
        self.F = np.array([[-1/self.G9,-1/self.G9],
                           [0,0]])
        self.L = np.array([[-1,-1,-1]])
        self.M = np.array([[1,0]])
        self.N = np.array([[0,0]])

        # Auxiliary Matrices
        self.Ix = np.eye(self.XDIM)
        self.Iv = np.eye(self.VDIM)
        self.H = np.linalg.inv(self.c*self.Ix - self.A)
        self.Ax = self.H @ (self.c*self.Ix + self.A)
        self.Bx = self.H @ self.B
        self.Cx = self.H @ self.C
        self.Ap = self.D @ self.H @ (self.c*self.Ix + self.A)
        self.Bp = self.D @ self.H @ self.B
        self.Cp = self.D @ self.H @ self.C
        self.K = self.D @ self.H @ self.C + self.F

        # DC Conditions
        self.x_dc = np.zeros((self.XDIM,1))
        self.u_dc = np.array([[self.Vcc]])
        self.v_dc = np.zeros((self.VDIM,1))
        self.computeInitialConditions()
        self.transistor.setOperatingPoint(self.v_dc)
        self.transistor.setLinearization(self.linearize)

        # Memory Vectors
        self.x_1 = self.x_dc.copy()
        self.u_1 = self.u.copy()
        self.v_1 = self.v_dc.copy()
        self.fv_1 = self.nonlinearFunction(self.v_1)
    
    # Public Methods
    def processSample(self,Vin):

        # Set Input Voltage
        self.u[0] = Vin

        # Solve Nonlinearities
        self.p = self.Ap @ self.x_1 + self.Bp @ (self.u + self.u_1) + self.Cp @ self.fv_1 + self.E @ self.u
        self.v = self.solveNonlinearity(self.v_1)
        self.fv = self.nonlinearFunction(self.v)

        # Update States
        self.x = self.Ax @ self.x_1 + self.Bx @ (self.u + self.u_1) + self.Cx @ (self.fv + self.fv_1)

        # Compute Output Voltage
        y = self.L @ self.x + self.M @ self.u + self.N @ self.fv
        Vout = y.item()

        # Update Memory
        self.x_1 = self.x
        self.u_1[0] = Vin
        self.v_1 = self.v
        self.fv_1 = self.fv

        # Return
        return Vout
    
    # Private Methods
    def computeInitialConditions(self):
        
        # DC System Matrices
        Req = self.R6 + self.R7 + self.R8
        A = np.array([[-self.R6/Req,-self.R7/Req,(self.R6 + self.R7)/Req]]).T
        B = np.array([[(self.R6*self.R7 + self.R6*self.R8)/Req,self.R6*self.R8/Req],
                      [-self.R6*self.R7/Req,self.R7*self.R8/Req],
                      [-self.R6*self.R8/Req,-(self.R8*self.R6 + self.R8*self.R7)/Req]])
        C = np.array([[self.R6/Req,-self.R7/Req]]).T
        D = np.array([[-(self.R6*self.R7 + self.R6*self.R8 + self.R9*self.R6 + self.R9*self.R7 + self.R9*self.R8)/Req,-(self.R6*self.R8 + self.R9*self.R6 + self.R9*self.R7 + self.R9*self.R8)/Req],
                      [-self.R6*self.R7/Req,self.R7*self.R8/Req]])
        
        # Define Auxiliary Functions
        residual_function = lambda v: C @ self.u_dc + D @ self.nonlinearFunction(v) - v
        residual_jacobian = lambda v: D @ self.nonlinearJacobian(v) - self.Iv

        # Compute Initial Conditions
        self.v_dc = nl.newtonRaphson(self.v_dc,residual_function,residual_jacobian,self.threshold,self.iterations,self.subiterations)
        self.x_dc = A @ self.u_dc + B @ self.nonlinearFunction(self.v_dc)
    
    def nonlinearFunction(self,v):

        # Compute Function Value
        f = self.transistor.current(v)

        # Return
        return f

    def nonlinearJacobian(self,v):

        # Compute Jacobian Value
        Jf = self.transistor.slope(v)
        
        # Return
        return Jf

    def residualFunction(self,v):

        # Compute Function Value
        g = self.p + self.K @ self.nonlinearFunction(v) - v

        # Return
        return g
    
    def residualJacobian(self,v):

        # Compute Jacobian Value
        Jg = self.K @ self.nonlinearJacobian(v) - self.Iv

        # Return
        return Jg
    
    def solveNonlinearity(self,v0):
        
        # Solve Nonlinearity
        v = nl.newtonRaphson(v0,self.residualFunction,self.residualJacobian,self.threshold,self.iterations,self.subiterations)

        # Return
        return v
    
# ====================================================================================================
# Class: OperationalAmplifier
# Description: This class represents the Boss DS-1 Operational Amplifier stage
# ====================================================================================================
class OperationalAmplifier(BaseModel):
    
    # Static Variables
    name = 'Operational Amplifier'
    label = 'WDF'
    parameters = ['distortion']

    # Constructor Method
    def __init__(self,fs = 44100,distortion = 0.5,**kwargs):

        # Initialize Base Class
        super().__init__(fs)

        # Resistors
        self.R11 = 100e3
        self.R13 = 4.7e3
        self.VR1 = 100e3
        self.distortion = utils.limit(distortion,0,1,utils.POT_MARGIN)
        self.Ra = self.VR1*(1 - self.distortion)
        self.Rb = self.VR1*self.distortion

        # Capacitors
        self.C7 = 100e-12
        self.C8 = 0.47e-6

        # Op-Amp
        self.Vsat = 4.5
        self.linearize = kwargs.get('linearize',False)

        # WDF Tree #1
        self.Vin_port = wdf.ResistiveVoltageSource(self.Ra)
        self.R13_port = wdf.Resistor(self.R13)
        self.C8_port = wdf.Capacitor(self.C8,self.fs)
        self.series_root = wdf.SeriesAdaptorRoot(self.Vin_port,self.R13_port,self.C8_port)
        
        # WDF Tree #2
        self.Ifb_port = wdf.ResistiveCurrentSource(self.Rb)
        self.C7_port = wdf.Capacitor(self.C7,self.fs)
        self.parallel_root = wdf.ParallelAdaptorRoot(self.Ifb_port,self.C7_port)

    # Public Methods
    def processSample(self,Vin):

        # Compute Feedback Current
        self.Vin_port.setVoltage(Vin)
        self.series_root.collectIncidentWaves()
        self.series_root.propagateReflectedWaves()
        Ifb = -self.Vin_port.waveToCurrent()

        # Compute Feedback Voltage
        self.Ifb_port.setCurrent(Ifb)
        self.parallel_root.collectIncidentWaves()
        self.parallel_root.propagateReflectedWaves()
        Vfb = self.Ifb_port.waveToVoltage()

        # Compute Output Voltage
        if self.linearize:
            Vout = Vfb + Vin
        else:
            Vout = utils.limit(Vfb + Vin,-self.Vsat,self.Vsat)

        # Return
        return Vout
    
    def setDistortion(self,distortion):

        # Check Change
        if(self.distortion != distortion):

            # Update Distortion Value
            self.distortion = utils.limit(distortion,0,1,utils.POT_MARGIN)

            # Update WDF Structure
            self.updateStructure()

    def getDistortion(self):
        return self.distortion

    # Private Methods
    def updateStructure(self):

        # Update Variable Resistors
        self.Ra = self.VR1*(1 - self.distortion)
        self.Rb = self.VR1*self.distortion

        # Update Port Resistances
        self.Vin_port.setPortResistance(self.Ra)
        self.Ifb_port.setPortResistance(self.Rb)

        # Update Scattering Matrices
        self.series_root.computeScatteringMatrix()
        self.parallel_root.computeScatteringMatrix()

# ====================================================================================================
# Class: DiodeClipper
# Description: This class represents the Boss DS-1 Diode Clipper stage
# ====================================================================================================
class DiodeClipper(BaseModel):
    
    # Static Variables
    name = 'Diode Clipper'
    label = 'PHS'
    parameters = []

    # Constructor Method
    def __init__(self,fs = 44100,**kwargs):

        # Initialize Base Class
        super().__init__(fs)

        # Dimension Variables
        self.SDIM = 2
        self.DDIM = 2
        self.PDIM = 3 - 1
        self.NDIM = 4 + 1
        self.BDIM = self.SDIM + self.DDIM + self.PDIM + 1

        # Resistors
        self.R14 = 2.2e3

        # Conductances
        self.G14 = 1/self.R14

        # Capacitors
        self.C9 = 0.47e-6
        self.C10 = 0.01e-6

        # Diodes
        self.Is = 2.52e-9
        self.n = 1.752
        self.Vt = utils.THERMAL_VOLTAGE
        self.diodes = nl.Diodes(Is = self.Is,n = self.n,Vt = self.Vt)
        self.linearize = kwargs.get('linearize',False)
        self.diodes.setLinearization(self.linearize)

        # Newton-Raphson Parameters
        self.threshold = utils.THRESHOLD
        self.iterations = utils.ITERATIONS
        self.subiterations = utils.SUBITERATIONS

        # System Vectors
        self.x = np.zeros((self.SDIM,1))
        self.w = np.zeros((self.DDIM,1))
        self.u = np.zeros((self.PDIM,1))
        self.p = np.zeros((self.DDIM,1))
        self.permutation = np.array([3,4,2,0,5,1,6])

        # System Matrices
        self.incident_matrix = np.array([[0,-1,0,0,0,-1,-1],
                                         [0,0,1,0,0,1,0],
                                         [0,0,-1,1,0,0,0],
                                         [1,1,0,-1,1,0,0],
                                         [-1,0,0,0,-1,0,1]])
        self.J = None
        self.computeMatrixJ()
        self.Jx = self.J[0:self.SDIM,0:self.SDIM]
        self.Jw = self.J[self.SDIM:(self.SDIM + self.DDIM),self.SDIM:(self.SDIM + self.DDIM)]
        self.Jy = self.J[(self.SDIM + self.DDIM):(self.SDIM + self.DDIM + self.PDIM),(self.SDIM + self.DDIM):(self.SDIM + self.DDIM + self.PDIM)]
        self.K = -self.J[0:self.SDIM,self.SDIM:(self.SDIM + self.DDIM)]
        self.Gx = -self.J[0:self.SDIM,(self.SDIM + self.DDIM):(self.SDIM + self.DDIM + self.PDIM)]
        self.Gw = -self.J[self.SDIM:(self.SDIM + self.DDIM),(self.SDIM + self.DDIM):(self.SDIM + self.DDIM + self.PDIM)]

        # Auxiliary Matrices
        self.Ix = np.eye(self.SDIM)
        self.Iw = np.eye(self.DDIM)
        self.Q = np.diag(np.array([1/self.C9,1/self.C10]))
        self.D = np.linalg.inv(self.Ix/self.Ts - (self.Jx @ self.Q)/2)
        self.Ax = self.D @ self.Jx @ self.Q
        self.Bx = -self.D @ self.K
        self.Cx = -self.D @ self.Gx
        self.Aw = (np.transpose(self.K) @ self.Q @ (2*self.Ix + self.Ax))/2
        self.Bw = self.Jw + (np.transpose(self.K) @ self.Q @ self.Bx)/2
        self.Cw = -self.Gw + (np.transpose(self.K) @ self.Q @ self.Cx)/2
        self.Ay = -((np.transpose(self.Gx) @ self.Q @ (2*self.Ix + self.Ax))/2)
        self.By = -(np.transpose(self.Gw) + (np.transpose(self.Gx) @ self.Q @ self.Bx)/2)
        self.Cy = -(self.Jy + (np.transpose(self.Gx) @ self.Q @ self.Cx)/2)

    # Public Methods
    def processSample(self,Vin):

        # Set Input Voltage
        self.u[0] = Vin

        # Solve Nonlinearities
        self.p = self.Aw @ self.x + self.Cw @ self.u
        self.w = self.solveNonlinearity(self.w)
        z = self.nonlinearFunction(self.w)

        # Compute Output Voltage
        y = self.Ay @ self.x + self.By @ z + self.Cy @ self.u
        Vout = -y[-1,0]

        # Update States
        dx = self.Ax @ self.x + self.Bx @ z + self.Cx @ self.u
        self.x += dx
        
        # Return
        return Vout
    
    # Private Methods
    def computeMatrixJ(self):

        # Compute Auxiliary Matrices
        P = np.eye(self.BDIM)[self.permutation]
        g1 = self.incident_matrix[1:,0:(self.BDIM - self.NDIM + 1)]
        g2 = self.incident_matrix[1:,(self.BDIM - self.NDIM + 1):]
        g = np.linalg.inv(g2) @ g1

        # Compute Matrix J
        J_tmp = np.zeros((self.BDIM,self.BDIM))
        J_tmp[(self.BDIM - g.shape[0]):,:g.shape[1]] = -g
        J_tmp[:g.shape[1],(self.BDIM - g.shape[0]):] = np.transpose(g)
        J_tmp = P @ J_tmp @ np.transpose(P)
        self.J = J_tmp[0:-1,0:-1]

    def nonlinearFunction(self,w):

        # Compute Function Value
        z = np.array([[w[0,0]*self.G14,self.diodes.current(w[1:]).item()]]).T

        # Return
        return z
    
    def nonlinearJacobian(self,w):

        # Compute Jacobian Value
        Jz = np.array([[self.G14,0],
                       [0,self.diodes.slope(w[1:]).item()]])
        
        # Return
        return Jz
    
    def residualFunction(self,w):

        # Compute Function Value
        g = self.p + self.Bw @ self.nonlinearFunction(w) - w

        # Return
        return g
    
    def residualJacobian(self,w):

        # Compute Jacobian Value
        Jg = self.Bw @ self.nonlinearJacobian(w) - self.Iw

        # Return
        return Jg
    
    def solveNonlinearity(self,w0):
        
        # Solve Nonlinearity
        w = nl.newtonRaphson(w0,self.residualFunction,self.residualJacobian,self.threshold,self.iterations,self.subiterations)

        # Return
        return w

# ====================================================================================================
# Class: PassiveFilter
# Description: This class represents the Boss DS-1 Passive Filter stage
# ====================================================================================================
class PassiveFilter(BaseModel):

    # Static Variables
    name = 'Passive Filter'
    label = 'IIR'
    parameters = ['tone']

    # Constructor Method
    def __init__(self,fs = 44100,tone = 0.5,**kwargs):

        # Initialize Base Class
        super().__init__(fs)

        # Resistors
        self.tone = utils.limit(tone,0,1)
        self.R15 = 2.2e3
        self.R16 = 6.8e3
        self.R17 = 6.8e3
        self.VR2 = 20e3
        self.VR3 = 100e3

        # Capacitors
        self.C11 = 22e-9
        self.C12 = 0.1e-6

        # System Memory
        self.order = 2
        self.x = np.zeros(self.order + 1)
        self.y = np.zeros(self.order)

        # System Coefficients
        self.exponents = np.arange(self.order + 1)
        self.k_B0 = self.k_B1 = self.k_B2 = self.b = []
        self.k_A0 = self.k_A1 = self.k_A2 = self.a = []
        self.initializeCoeffs()

    # Public Methods
    def processSample(self,Vin):

        # Set Input Voltage
        self.x[0] = Vin

        # Compute Output Voltage
        Vout = np.inner(self.b,self.x) - np.inner(self.a,self.y)

        # Update Memory
        for n in range(self.order,0,-1):
            self.x[n] = self.x[n - 1]
            if n > 1:
                self.y[n - 1] = self.y[n - 2]
        self.y[0] = Vout
        
        # Return
        return Vout
    
    def setTone(self,tone):

        # Check Change
        if(self.tone != tone):

            # Update Tone Value
            self.tone = utils.limit(tone,0,1)

            # Update Coefficients
            self.updateCoeffs()

    def getTone(self):
        return self.tone

    # Private Methods
    def initializeCoeffs(self):
        
        # Coefficient B0
        k_B0_0 = -self.C11*self.R15*self.R17*self.VR3*self.c - self.C11*self.R15*self.VR3*self.VR2*self.c - self.C11*self.R16*self.R17*self.VR3*self.c - self.C11*self.R17*self.VR3*self.VR2*self.c - self.R17*self.VR3 - self.VR3*self.VR2
        k_B0_1 = -self.C11*self.C12*self.R16*self.R17*self.VR3*self.VR2*np.power(self.c,2) + self.C11*self.R15*self.VR3*self.VR2*self.c + self.VR3*self.VR2
        k_B0_2 = 0
        self.k_B0 = np.array([k_B0_0,k_B0_1,k_B0_2])

        # Coefficient B1
        k_B1_0 = -2*self.R17*self.VR3 - 2*self.VR3*self.VR2
        k_B1_1 = 2*self.C11*self.C12*self.R16*self.R17*self.VR3*self.VR2*np.power(self.c,2) + 2*self.VR3*self.VR2
        k_B1_2 = 0
        self.k_B1 = np.array([k_B1_0,k_B1_1,k_B1_2])

        # Coefficient B2
        k_B2_0 = self.C11*self.R15*self.R17*self.VR3*self.c + self.C11*self.R15*self.VR3*self.VR2*self.c + self.C11*self.R16*self.R17*self.VR3*self.c + self.C11*self.R17*self.VR3*self.VR2*self.c - self.R17*self.VR3 - self.VR3*self.VR2
        k_B2_1 = -self.C11*self.C12*self.R16*self.R17*self.VR3*self.VR2*np.power(self.c,2) - self.C11*self.R15*self.VR3*self.VR2*self.c + self.VR3*self.VR2
        k_B2_2 = 0
        self.k_B2 = np.array([k_B2_0,k_B2_1,k_B2_2])

        # Coefficient Α0
        k_A0_0 = -self.C11*self.C12*self.R15*self.R16*self.R17*self.VR3*np.power(self.c,2) - self.C11*self.C12*self.R15*self.R16*self.VR3*self.VR2*np.power(self.c,2) - self.C11*self.C12*self.R16*self.R17*self.VR3*self.VR2*np.power(self.c,2) - self.C11*self.R15*self.R16*self.R17*self.c - self.C11*self.R15*self.R16*self.VR3*self.c - self.C11*self.R15*self.R16*self.VR2*self.c - self.C11*self.R15*self.R17*self.VR3*self.c - self.C11*self.R15*self.VR3*self.VR2*self.c - self.C11*self.R16*self.R17*self.VR3*self.c - self.C11*self.R16*self.R17*self.VR2*self.c - self.C11*self.R17*self.VR3*self.VR2*self.c - self.C12*self.R16*self.R17*self.VR3*self.c - self.C12*self.R16*self.VR3*self.VR2*self.c - self.R16*self.R17 - self.R16*self.VR3 - self.R16*self.VR2 - self.R17*self.VR3 - self.VR3*self.VR2
        k_A0_1 = -self.C11*self.C12*self.R15*self.R16*self.R17*self.VR2*np.power(self.c,2) - self.C11*self.C12*self.R15*self.R16*np.power(self.VR2,2)*np.power(self.c,2) - self.C11*self.C12*self.R16*self.R17*np.power(self.VR2,2)*np.power(self.c,2) + self.C11*self.R15*self.R16*self.VR2*self.c - self.C11*self.R15*self.R17*self.VR2*self.c - self.C11*self.R15*np.power(self.VR2,2)*self.c + self.C11*self.R16*self.R17*self.VR2*self.c - self.C11*self.R17*np.power(self.VR2,2)*self.c - self.C12*self.R16*self.R17*self.VR2*self.c - self.C12*self.R16*np.power(self.VR2,2)*self.c + self.R16*self.VR2 - self.R17*self.VR2 - np.power(self.VR2,2)
        k_A0_2 = self.C11*self.C12*self.R15*self.R16*np.power(self.VR2,2)*np.power(self.c,2) + self.C11*self.C12*self.R16*self.R17*np.power(self.VR2,2)*np.power(self.c,2) + self.C11*self.R15*np.power(self.VR2,2)*self.c + self.C11*self.R17*np.power(self.VR2,2)*self.c + self.C12*self.R16*np.power(self.VR2,2)*self.c + np.power(self.VR2,2)
        self.k_A0 = np.array([k_A0_0,k_A0_1,k_A0_2])
        
        # Coefficient Α1
        k_A1_0 = 2*self.C11*self.C12*self.R15*self.R16*self.R17*self.VR3*np.power(self.c,2) + 2*self.C11*self.C12*self.R15*self.R16*self.VR3*self.VR2*np.power(self.c,2) + 2*self.C11*self.C12*self.R16*self.R17*self.VR3*self.VR2*np.power(self.c,2) - 2*self.R16*self.R17 - 2*self.R16*self.VR3 - 2*self.R16*self.VR2 - 2*self.R17*self.VR3 - 2*self.VR3*self.VR2
        k_A1_1 = 2*self.C11*self.C12*self.R15*self.R16*self.R17*self.VR2*np.power(self.c,2) + 2*self.C11*self.C12*self.R15*self.R16*np.power(self.VR2,2)*np.power(self.c,2) + 2*self.C11*self.C12*self.R16*self.R17*np.power(self.VR2,2)*np.power(self.c,2) + 2*self.R16*self.VR2 - 2*self.R17*self.VR2 - 2*np.power(self.VR2,2)
        k_A1_2 = -2*self.C11*self.C12*self.R15*self.R16*np.power(self.VR2,2)*np.power(self.c,2) - 2*self.C11*self.C12*self.R16*self.R17*np.power(self.VR2,2)*np.power(self.c,2) + 2*np.power(self.VR2,2)
        self.k_A1 = np.array([k_A1_0,k_A1_1,k_A1_2])
        
        # Coefficient Α2
        k_A2_0 = -self.C11*self.C12*self.R15*self.R16*self.R17*self.VR3*np.power(self.c,2) - self.C11*self.C12*self.R15*self.R16*self.VR3*self.VR2*np.power(self.c,2) - self.C11*self.C12*self.R16*self.R17*self.VR3*self.VR2*np.power(self.c,2) + self.C11*self.R15*self.R16*self.R17*self.c + self.C11*self.R15*self.R16*self.VR3*self.c + self.C11*self.R15*self.R16*self.VR2*self.c + self.C11*self.R15*self.R17*self.VR3*self.c + self.C11*self.R15*self.VR3*self.VR2*self.c + self.C11*self.R16*self.R17*self.VR3*self.c + self.C11*self.R16*self.R17*self.VR2*self.c + self.C11*self.R17*self.VR3*self.VR2*self.c + self.C12*self.R16*self.R17*self.VR3*self.c + self.C12*self.R16*self.VR3*self.VR2*self.c - self.R16*self.R17 - self.R16*self.VR3 - self.R16*self.VR2 - self.R17*self.VR3 - self.VR3*self.VR2
        k_A2_1 = -self.C11*self.C12*self.R15*self.R16*self.R17*self.VR2*np.power(self.c,2) - self.C11*self.C12*self.R15*self.R16*np.power(self.VR2,2)*np.power(self.c,2) - self.C11*self.C12*self.R16*self.R17*np.power(self.VR2,2)*np.power(self.c,2) - self.C11*self.R15*self.R16*self.VR2*self.c + self.C11*self.R15*self.R17*self.VR2*self.c + self.C11*self.R15*np.power(self.VR2,2)*self.c - self.C11*self.R16*self.R17*self.VR2*self.c + self.C11*self.R17*np.power(self.VR2,2)*self.c + self.C12*self.R16*self.R17*self.VR2*self.c + self.C12*self.R16*np.power(self.VR2,2)*self.c + self.R16*self.VR2 - self.R17*self.VR2 - np.power(self.VR2,2)
        k_A2_2 = self.C11*self.C12*self.R15*self.R16*np.power(self.VR2,2)*np.power(self.c,2) + self.C11*self.C12*self.R16*self.R17*np.power(self.VR2,2)*np.power(self.c,2) - self.C11*self.R15*np.power(self.VR2,2)*self.c - self.C11*self.R17*np.power(self.VR2,2)*self.c - self.C12*self.R16*np.power(self.VR2,2)*self.c + np.power(self.VR2,2)
        self.k_A2 = np.array([k_A2_0,k_A2_1,k_A2_2])

        # Coefficient Initialization
        self.updateCoeffs()

    def updateCoeffs(self):

        # Coefficient Evaluation
        powers = self.tone**self.exponents
        B0 = np.inner(self.k_B0,powers)
        B1 = np.inner(self.k_B1,powers)
        B2 = np.inner(self.k_B2,powers)
        A0 = np.inner(self.k_A0,powers)
        A1 = np.inner(self.k_A1,powers)
        A2 = np.inner(self.k_A2,powers)

        # Coefficient Normalization
        self.b = np.array([B0,B1,B2])/A0
        self.a = np.array([A1,A2])/A0

# ====================================================================================================
# Class: BossDS1
# Description: This class represents the complete Boss DS-1 virtual analog model
# ====================================================================================================
class BossDS1(BaseModel):
    
    # Static Variables
    name = 'Boss DS-1'
    label = 'VAM'
    parameters = ['distortion','tone','level']

    # Constructor Method
    def __init__(self,fs = 44100,distortion = 0.5,tone = 0.5,level = 1,**kwargs):

        # Initialize Base Class
        super().__init__(fs)

        # Assign Arguments
        self.level = level
        self.linearize = kwargs.get('linearize',False)

        # Define Member Variables
        self.transistor_amplifier = TransistorAmplifier(fs = self.fs,linearize = self.linearize)
        self.operational_amplifier = OperationalAmplifier(fs = self.fs,distortion = distortion,linearize = self.linearize)
        self.diode_clipper = DiodeClipper(fs = self.fs,linearize = self.linearize)
        self.passive_filter = PassiveFilter(fs = self.fs,tone = tone)
        self.distortion = self.operational_amplifier.getDistortion()
        self.tone = self.passive_filter.getTone()

    # Public Methods
    def processSample(self,Vin):

        # Compute Output Voltage
        Vout = Vin
        Vout = self.transistor_amplifier.processSample(Vout)
        Vout = self.operational_amplifier.processSample(Vout)
        Vout = self.diode_clipper.processSample(Vout)
        Vout = self.passive_filter.processSample(Vout)
        Vout = self.level*Vout

        # Return
        return Vout
    
    def setDistortion(self,distortion):

        # Update Distortion Value
        self.operational_amplifier.setDistortion(distortion)

    def setTone(self,tone):

        # Update Tone Value
        self.passive_filter.setTone(tone)

    def setLevel(self,level):

        # Check Change
        if(self.level != level):

            # Update Level Value
            self.level = utils.limit(level,0,1)

    def getDistortion(self):
        return self.distortion
    
    def getTone(self):
        return self.tone
    
    def getLevel(self):
        return self.level