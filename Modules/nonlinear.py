# ====================================================================================================
# # ==================================================================================================
# # Filename: nonlinear.py
# #
# # Description: This module contains functions and classes related to nonlinear circuits
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# Import Modules (Built-In)
from abc import ABC,abstractmethod

# Import Modules (User-Defined)
import utilities as utils

# Import Modules (Third-Party)
import numpy as np

# ====================================================================================================
# Function: newtonRaphson()
# Description: This function iteratively solves a nonlinear equation via the Newton-Raphson method
# ====================================================================================================
def newtonRaphson(x0,function,jacobian,threshold = utils.THRESHOLD,iterations = utils.ITERATIONS,subiterations = utils.SUBITERATIONS):
    
    # Initialize Variables
    iter = 0
    x = x0
    g = function(x)
    norm = np.inner(g.flatten(),g.flatten())
    error = np.max(np.abs(g))

    # Perform Iterative Process
    while(iter < iterations and error > threshold):

        # Initialize Iteration
        J = jacobian(x)
        step = np.linalg.solve(J,g)
        x -= step
        g = function(x)
        norm_tmp = np.inner(g.flatten(),g.flatten())

        # Apply Damping
        subiter = 0
        while(subiter < subiterations and (norm_tmp > norm or not np.isfinite(norm_tmp))):
            step /= 2
            x += step
            g = function(x)
            norm_tmp = np.inner(g.flatten(),g.flatten())
            subiter += 1

        # Conclude Iteration
        norm = norm_tmp
        error = np.max(np.abs(g))
        iter += 1

    # Return
    return x

# ====================================================================================================
# Class: BaseNonlinear
# Description: This class serves as a base for nonlinear component classes
# ====================================================================================================
class BaseNonlinear(ABC):

    # Constructor Method
    @abstractmethod
    def __init__(self,op = None,linearize = False):

        # Assign Arguments
        self.op = op
        self.linearize = linearize
    
    # Public Methods
    @abstractmethod
    def current(self,v):
        pass

    @abstractmethod
    def slope(self,v):
        pass

    def setOperatingPoint(self,op):
        self.op = op

    def setLinearization(self,linearize):
        self.linearize = linearize

# ====================================================================================================
# Class: Diodes
# Description: This class represents a pair of anti-parallel diodes
# ====================================================================================================
class Diodes(BaseNonlinear):

    # Constructor Method
    def __init__(self,Is,n = 1,Vt = utils.THERMAL_VOLTAGE,op = np.zeros((1,1)),linearize = False):

        # Initialize Base Class
        super().__init__(op,linearize)

        # Assign Arguments
        self.Is = Is
        self.n = n
        self.Vt = Vt

        # Define Member Variables
        self.k = 1/(self.n*self.Vt)
        self.d = 2*self.Is
        self.kd = self.k*self.d

    # Public Method
    def current(self,v):

        # Calculate Current
        if self.linearize:
            i = self.d*np.sinh(self.op[0,0]*self.k) + self.kd*np.cosh(self.op[0,0]*self.k)*(v[0,0] - self.op[0,0])
        else:
            i = self.d*np.sinh(v[0,0]*self.k)
        current = np.array([[i]])

        # Return
        return current
    
    def slope(self,v):

        # Calculate Slope
        if self.linearize:
            s = self.kd*np.cosh(self.op[0,0]*self.k)
        else:
            s = self.kd*np.cosh(v[0,0]*self.k)
        slope = np.array([[s]])
        
        # Return
        return slope

# ====================================================================================================
# Class: Transistor
# Description: This class represents a BJT transistor component
# ====================================================================================================
class Transistor(BaseNonlinear):

    # Constructor Method
    def __init__(self,Is,Bf,Br,n = 1,Vt = utils.THERMAL_VOLTAGE,op = np.zeros((2,1)),linearize = False):
        
        # Initialize Base Class
        super().__init__(op,linearize)

        # Assign Arguments
        self.Bf = Bf
        self.Br = Br
        self.Is = Is
        self.n = n
        self.Vt = Vt

        # Define Member Variables
        self.k = 1/(self.n*self.Vt)
        self.b0 = self.Is/self.Bf
        self.b1 = self.Is/self.Br
        self.c0 = self.Is
        self.c1 = -self.Is/(self.Br/(1 + self.Br))
        self.kb0 = self.k*self.b0
        self.kb1 = self.k*self.b1
        self.kc0 = self.k*self.c0
        self.kc1 = self.k*self.c1

    # Public Method
    def current(self,v):

        # Calculate Current
        if self.linearize:
            exp0 = np.exp(self.op[0,0]*self.k) - 1
            exp1 = np.exp(self.op[1,0]*self.k) - 1
            exp00 = np.exp(self.op[0,0]*self.k)*(v[0,0] - self.op[0,0])
            exp10 = np.exp(self.op[1,0]*self.k)*(v[1,0] - self.op[1,0])
            ib = (self.b0*exp0 + self.b1*exp1) + self.kb0*exp00 + self.kb1*exp10
            ic = (self.c0*exp0 + self.c1*exp1) + self.kc0*exp00 + self.kc1*exp10
        else:
            exp0 = np.exp(v[0,0]*self.k) - 1
            exp1 = np.exp(v[1,0]*self.k) - 1
            ib = self.b0*exp0 + self.b1*exp1
            ic = self.c0*exp0 + self.c1*exp1
        current = np.array([[ib,ic]]).T

        # Return
        return current
    
    def slope(self,v):

        # Calculate Slope
        if self.linearize:
            exp0 = np.exp(self.op[0,0]*self.k)
            exp1 = np.exp(self.op[1,0]*self.k)
            sb0 = self.kb0*exp0
            sb1 = self.kb1*exp1
            sc0 = self.kc0*exp0
            sc1 = self.kc1*exp1
        else:
            exp0 = np.exp(v[0,0]*self.k)
            exp1 = np.exp(v[1,0]*self.k)
            sb0 = self.kb0*exp0
            sb1 = self.kb1*exp1
            sc0 = self.kc0*exp0
            sc1 = self.kc1*exp1
        slope = np.array([[sb0,sb1],
                          [sc0,sc1]])
        
        # Return
        return slope