# SPDX-FileCopyrightText: 2024 Bernd Flemisch, Dennis Gläser <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

""" Generate table for CO2 - H2O solubilities.

The calculations implement formulas (11)-(14) from
Spycher, Pruess, Ennis-King (2003),
https://doi.org/10.1016/S0016-7037(03)00273-4.
"""

import argparse
import math
import numpy as np  # pylint: disable=import-error
import matplotlib.pyplot as plt

try:
    import tqdm
    def withProgress(iterable):
        return tqdm.tqdm(iterable)
except ImportError:
    def withProgress(iterable):
        return iterable


class ParameterRange:
    def __init__(self, min: float, max: float, numSamples: int) -> None:
        assert (max > min and numSamples > 1) or (max == min and numSamples == 1)
        self._min = min
        self._max = max
        self._numSamples = numSamples

    def __getitem__(self, i: int) -> float:
        return self._min + i*self.step

    @property
    def min(self) -> float:
        return self._min

    @property
    def max(self) -> float:
        return self._max

    @property
    def numSamples(self) -> float:
        return self._numSamples

    @property
    def step(self) -> float:
        if self._numSamples == 1:
            return 0
        else:
            return (self._max - self._min)/float(self._numSamples - 1)

def equilibriumConstantCO2(T):
    TinC = T - 273.15 # temperature in °C
    c = [1.189, 1.304e-2, -5.446e-5]
    logk0_CO2 = c[0] + c[1] * TinC + c[2] * TinC * TinC
    return math.pow(10, logk0_CO2)

def equilibriumConstantH2O(T):
    TinC = T - 273.15 # temperature in °C
    c = [-2.209, 3.097e-2, -1.098e-4, 2.048e-7]
    logk0_H2O = c[0] + c[1]*TinC + c[2]*TinC*TinC + c[3]*TinC*TinC*TinC
    return math.pow(10, logk0_H2O)

# Adapted from https://github.com/mwburgoyne/pyResToolbox
def molarVolume(T, p):
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]
    p_bar = p/1e5 # phase pressure in bar
    RTp = R*T/p_bar
    a_CO2 = (7.54e7 - 4.13e4*T) # mixture parameter of Redlich-Kwong equation
    b_CO2 = 27.8 # mixture parameter of Redlich-Kwong equation
    aT12p = a_CO2/(p_bar*T**0.5)
    
    # calculate coefficients for cubic equation
    e2 = -RTp
    e1 = -(RTp*b_CO2 - aT12p + b_CO2**2)
    e0 = -aT12p*b_CO2

    # solve the cubic equation    
    Z = np.roots(np.array([1.0, e2, e1, e0]))
    Z = np.array([x for x in Z if np.isreal(x)]) # Keep only real results
    if len(Z) > 1: # Evaluate which root to use per Eqs 25 and 26 in Spycher & Pruess (2003)
        vgas, vliq = max(Z), min(Z)
            
        w1 = p_bar*(vgas - vliq)
        w2 = R*T*np.log((vgas - b_CO2)/(vliq - b_CO2)) + a_CO2/(T**0.5*b_CO2)*np.log((vgas + b_CO2)*vliq/((vliq + b_CO2)*vgas))
            
        if w2 - w1 > 0:
            Z[0] = max(Z)
        else:
            Z[0] = min(Z)
                
    return np.real(Z[0])

def fugacityCoefficientCO2(T, p):
    V = molarVolume(T, p)
    p_bar = p/1e5 # phase pressure in bar
    a_CO2 = (7.54e7 - 4.13e4*T) # mixture parameter of Redlich-Kwong equation
    b_CO2 = 27.8 # mixture parameter of Redlich-Kwong equation
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]

    lnPhiCO2 = math.log(V/(V - b_CO2)) + b_CO2/(V - b_CO2) \
               - 2*a_CO2/(R*math.pow(T, 1.5)*b_CO2)*math.log((V + b_CO2)/V) \
               + a_CO2*b_CO2/(R*math.pow(T, 1.5)*b_CO2*b_CO2) \
                 *(math.log((V + b_CO2)/V) - b_CO2/(V + b_CO2)) \
               - math.log(p_bar*V/(R*T))
    return math.exp(lnPhiCO2)

def fugacityCoefficientH2O(T, p):
    V = molarVolume(T, p)
    p_bar = p/1e5 # phase pressure in bar
    a_CO2 = (7.54e7 - 4.13e4*T) # mixture parameter of  Redlich-Kwong equation
    a_CO2_H2O = 7.89e7 # mixture parameter of Redlich-Kwong equation
    b_CO2 = 27.8 # mixture parameter of Redlich-Kwong equation
    b_H2O = 18.18 # mixture parameter of Redlich-Kwong equation
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]

    lnPhiH2O = math.log(V/(V - b_CO2)) + b_H2O/(V - b_CO2) \
               - 2*a_CO2_H2O/(R*math.pow(T, 1.5)*b_CO2)*math.log((V + b_CO2)/V) \
               + a_CO2*b_H2O/(R*math.pow(T, 1.5)*b_CO2*b_CO2) \
                 *(math.log((V + b_CO2)/V) - b_CO2/(V + b_CO2)) \
               - math.log(p_bar*V/(R*T))
    return math.exp(lnPhiH2O)

def computeA(T, p):
    deltaP = p/1e5 - 1 # pressure range [bar] from p0 = 1 bar to p
    v_av_H2O = 18.1 # average partial molar volume of H2O [cm3/mol]
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]
    k0_H2O = equilibriumConstantH2O(T) # equilibrium constant for H2O at 1 bar
    phi_H2O = fugacityCoefficientH2O(T, p) # fugacity coefficient of H2O for the water-CO2 system
    p_bar = p/1e5
    return k0_H2O/(phi_H2O*p_bar)*math.exp(deltaP*v_av_H2O/(R*T))

def computeB(T, p):
    deltaP = p/1e5 - 1 # pressure range [bar] from p0 = 1 bar to p
    v_av_CO2 = 32.6 # average partial molar volume of CO2 [cm3/mol]
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]
    k0_CO2 = equilibriumConstantCO2(T) # equilibrium constant for CO2 at 1 bar
    phi_CO2 = fugacityCoefficientCO2(T, p) # fugacity coefficient of CO2 for the water-CO2 system
    p_bar = p/1e5
    return phi_CO2*p_bar/(55.508*k0_CO2)*math.exp(-(deltaP*v_av_CO2)/(R*T))

def makeSolubilityTable():
    parser = argparse.ArgumentParser(
        description="This script generates tables for CO2 - H2O solubilities \n"
        "according to Spycher, Pruess, Ennis-King (2003).\n"
    )
    parser.add_argument(
        "-t1", "--min_temp", required=True, type=float, help="The minimum temperature in degree Celcius."
    )
    parser.add_argument(
        "-t2", "--max_temp", required=True, type=float, help="The maximum temperature in degree Celcius."
    )
    parser.add_argument(
        "-nt",
        "--n_temp",
        required=True,
        type=int,
        help="The number of temperature sampling points."
        " min_temp is the first sampling point, max_temp the last.",
    )
    parser.add_argument(
        "-p1", "--min_press", required=True, type=float, help="The minimum phase pressure in Pascal."
    )
    parser.add_argument(
        "-p2", "--max_press", required=True, type=float, help="The maximum phase pressure in Pascal."
    )
    parser.add_argument(
        "-np",
        "--n_press",
        required=True,
        type=int,
        help="The number of pressure sampling points."
        " min_press is the first sampling point, max_press the last.",
    )
    parser.add_argument(
        "--visualize", required=False, help="Set to true for visualizing the results.", action=argparse.BooleanOptionalAction
    )
    parser.set_defaults(visualize=False)

    cmdArgs = vars(parser.parse_args())

    vis = cmdArgs["visualize"]
    pressures = ParameterRange(
        min=cmdArgs["min_press"],
        max=cmdArgs["max_press"],
        numSamples=cmdArgs["n_press"]
    )
    temperatures = ParameterRange(
        min=cmdArgs["min_temp"],
        max=cmdArgs["max_temp"],
        numSamples=cmdArgs["n_temp"]
    )

    fileName = "solubilities.csv"
    outFile = open(fileName, "w")
    outFile.write(f"# This autogenerated file contains solubilities of CO2 and H2O in a respective fluid system.\n")
    outFile.write("# The values have been calculated by means of (11)-(14) in https://doi.org/10.1016/S0016-7037(03)00273-4.\n#\n")
    outFile.write("# Concerning temperature and pressure ranges, the following parameters have been used:\n")
    outFile.write(f"# min temperature = {temperatures.min}, max temperature = {temperatures.max}, #temperature sampling points = {temperatures.numSamples}\n")
    outFile.write(f"# min phase pressure = {pressures.min}, max phase pressure = {pressures.max}, #pressure sampling points = {pressures.numSamples}\n#\n")
    outFile.write("# temperature [°C], phase pressure [Pa],   y_H2O [mol/mol],   x_CO2 [mol/mol]\n")

    # In the following, the notation and equation numbers
    # from Spycher et al. 2003 are used.
    for i in range(temperatures.numSamples):
        T = temperatures[i]
        yVec = []
        xVec = []
        pVec = []
        for j in range(pressures.numSamples):
            p = pressures[j]
            pVec.append(p/1e5)
            # convert temperature to Kelvin for the function calls
            A = computeA(T + 273.15, p); # (11)
            B = computeB(T + 273.15, p); # (12)
            y_H2O = (1 - B)/(1/A - B) # (13)
            x_CO2 = B*(1 - y_H2O) # (14)
            yVec.append(1e3*y_H2O)
            xVec.append(1e2*x_CO2)

            outFile.write(f" {T:.11e},   {p:.11e}, {y_H2O:.11e}, {x_CO2:.11e}\n")

        if vis:
            ax = plt.subplot(2, 1, 1)
            ax.plot(pVec, yVec)
            ax.set_title(f'y_H2O at {T} °C')
            ax.set_xlabel(f'p [bar]')
            ax.set_ylabel(f'y_H2O x 1000')
            ax.set_ylim((0, yVec[-1]+3))
            plt.grid(visible=True)
            ax = plt.subplot(2, 1, 2)
            ax.plot(pVec, xVec)
            ax.set_title(f'x_CO2 at {T} °C')
            ax.set_xlabel(f'p [bar]')
            ax.set_ylabel(f'x_CO2 x 100')
            ax.set_ylim((0, xVec[-1]+2))
            plt.grid(visible=True)
            plt.show()

    print(f"A file {fileName} has been generated.")

if __name__ == "__main__":
    makeSolubilityTable()
