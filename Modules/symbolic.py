# ====================================================================================================
# # ==================================================================================================
# # Filename: symbolic.py
# #
# # Description: This module contains functions used for symbolic computations
# #
# # Author: Alexandros Iliadis
# # Project: MSc Thesis
# # Date: June 2025
# # ==================================================================================================
# ====================================================================================================

# Import Modules (Built-In)
from collections.abc import Iterable

# Import Modules (User-Defined)
None

# Import Modules (Third-Party)
import sympy as sp

# ====================================================================================================
# Function: systemToMatrixForm()
# Description: This function converts a system of linear equations into matrix form
# ====================================================================================================
def systemToMatrixForm(equations,x):

    # Derive System Matrices
    A,b = sp.linear_eq_to_matrix(equations,x)

    # Adjust Signs
    for row in range(len(b)):
        if b[row].could_extract_minus_sign():
            b[row] *= -1
            A[row,:] *= -1

    # Return
    return A,b

# ====================================================================================================
# Function: solveLinearSystem()
# Description: This function solves a linear system given in matrix form
# ====================================================================================================
def solveLinearSystem(A,b,x):

    # Solve Linear System
    solution = sp.linsolve((A,b),x).args[0]
    solution = list(sp.together(solution))

    # Return
    return solution

# ====================================================================================================
# Function: formEquations()
# Description: This function combines LHS and RHS expressions to form a list of equations
# ====================================================================================================
def formEquations(lhs,rhs):

    # Check Argument Types
    if not isinstance(lhs,Iterable):
        lhs = [lhs]
    if not isinstance(rhs,Iterable):
        rhs = [rhs]

    # Check Argument Lengths
    if len(lhs) != len(rhs):
        raise ValueError("LHS and RHS lists must have the same length")
    
    # Form List of Equations
    equations = [sp.Eq(l,r) for l,r in zip(lhs,rhs)]

    # Return
    return equations

# ====================================================================================================
# Function: splitEquations()
# Description: This function splits a list of equations into two lists of LHS and RHS expressions
# ====================================================================================================
def splitEquations(equations):

    # Check Argument Type
    if not isinstance(equations,Iterable):
        equations = [equations]

    # Initialize Lists
    lhs = []
    rhs = []

    # Split LHS & RHS
    for eq in equations:
        lhs.append(eq.lhs)
        rhs.append(eq.rhs)

    # Return
    return lhs,rhs

# ====================================================================================================
# Function: extractEquations()
# Description: This function extracts a sublist from a list of equations based on specified LHS symbols
# ====================================================================================================
def extractEquations(equations,symbols):

    # Check Argument Type
    if not isinstance(symbols,Iterable):
        symbols = [symbols]

    # Initialize Sublist
    sublist = []

    # Extract Sublist
    for eq in equations:
        if eq.lhs in symbols:
            sublist.append(eq)

    # Return
    return sublist

# ====================================================================================================
# Function: substituteSymbols()
# Description: This function replaces symbols in an expression based on their respective equations
# ====================================================================================================
def substituteSymbols(expression,equations):

    # Check Argument Type
    if not isinstance(equations,Iterable):
        equations = [equations]
    
    # Replace Symbols
    expression = expression.subs([(eq.lhs,eq.rhs) for eq in equations])

    # Return
    return expression

# ====================================================================================================
# Function: polyCoeffs()
# Description: This function retrieves the coefficients from a polynomial expression
# ====================================================================================================
def polyCoeffs(polynomial,coeff_prefix = 'c'):

    # Create Coefficient Symbols
    degree = sp.degree(polynomial)
    lhs = sp.symbols(f'{coeff_prefix}_0:{degree + 1}')

    # Retrieve Polynomial Coefficients
    rhs = list(reversed(polynomial.all_coeffs()))

    # Form Coefficient Equations
    coeff_eq = formEquations(lhs,rhs)

    # Return
    return coeff_eq

# ====================================================================================================
# Function: rationalFraction()
# Description: This function arranges an expression into a rational fraction with explicitly defined coefficients
# ====================================================================================================
def rationalFraction(expression,symbol,num_prefix = 'b',den_prefix = 'a'):

    # Split Rational Fraction
    num,den = sp.fraction(sp.cancel(sp.together(expression)))

    # Extract Numerator Coefficients
    num_poly = sp.Poly(num,symbol)
    num_coeffs_eq = polyCoeffs(num_poly,num_prefix)

    # Extract Denominator Coefficients
    den_poly = sp.Poly(den,symbol)
    den_coeffs_eq = polyCoeffs(den_poly,den_prefix)

    # Form Rational Fraction
    rf_num = sum(num_coeffs_n*symbol**n for n,num_coeffs_n in enumerate(splitEquations(num_coeffs_eq)[0]))
    rf_den = sum(den_coeffs_n*symbol**n for n,den_coeffs_n in enumerate(splitEquations(den_coeffs_eq)[0]))
    rf = rf_num/rf_den

    # Return
    return rf,num_coeffs_eq,den_coeffs_eq

# ====================================================================================================
# Function: groupBySymbol()
# Description: This function groups the terms of an expression with respect to given a symbol
# ====================================================================================================
def groupBySymbol(expression,symbol):

    # Group Expression Terms
    expression = sp.collect(sp.expand(expression),symbol)

    # Return
    return expression

# ====================================================================================================
# Function: kMethodMatrices()
# Description: This function derives the K-Method matrices for a given list of state-space equations
# ====================================================================================================
def kMethodMatrices(equations,x,u,i):

    # Initialize Matrices
    A = sp.zeros(len(equations),len(x))
    B = sp.zeros(len(equations),len(u))
    C = sp.zeros(len(equations),len(i))

    # Fill Matrices
    for row in range(len(equations)):
        expression = sp.expand(equations[row].rhs)

        # Matrix A
        for col in range(len(x)):
            A[row,col] = expression.coeff(x[col])

        # Matrix B
        for col in range(len(u)):
            B[row,col] = expression.coeff(u[col]) 

        # Matrix C
        for col in range(len(i)):
            C[row,col] = expression.coeff(i[col]) 

    A = sp.together(A)
    B = sp.together(B)
    C = sp.together(C)
    
    # Return
    return A,B,C

# ====================================================================================================
# Function: kMethodMatricesDC()
# Description: This function derives the K-Method matrices for a given list of DC state-space equations
# ====================================================================================================
def kMethodMatricesDC(equations,u,i):

    # Initialize Matrices
    A = sp.zeros(len(equations),len(u))
    B = sp.zeros(len(equations),len(i))

    # Fill Matrices
    for row in range(len(equations)):
        expression = sp.expand(equations[row].rhs)

        # Matrix A
        for col in range(len(u)):
            A[row,col] = expression.coeff(u[col]) 

        # Matrix B
        for col in range(len(i)):
            B[row,col] = expression.coeff(i[col])

    A = sp.together(A)
    B = sp.together(B)
    
    # Return
    return A,B