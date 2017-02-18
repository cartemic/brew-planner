# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 10:37:56 2017

@author: Mick
"""

# %%  BREW PLANNER
# under construction
# for planning and documenting homebrews


# %% IMPORT PACKAGES
import numpy as np
import matplotlib.pyplot as plt


# %% DEFINE TERMS
# GOAL: build a GUI for these inputs


# %% DEFINE FUNCTIONS


def correctSG(SG_reading, T):
    # Corrects specific gravity readings for temperature
    #
    # INPUTS:
    # SG_reading:   SG as read from hydrometer
    # T:            Fluid temperature at time of reading (F)
    #
    # OUTPUTS:
    # SG:           SG corrected for temperature
    
    # SG correction factor coefficients
    a = 0.00130346
    b = 1.34722124e-4
    c = 2.04052596e-6
    d = 2.32820948e-9
    
    # Calculate SG correction factor from third order polyfit
    SG_corr = a - b * T + c * T**2 - d * T**3
    
    # Correct SG reading
    SG = SG_reading + SG_corr
    
    # Output corrected reading
    return(SG)


def predictColor(grainBill, batchVol):
    # Predicts final color (SRM) using the Morey formula
    #
    # INPUTS:
    # grainBill:    dictionary containing at least
    #   L:  array of grain colors (Lovibond)
    #   m:  array of grain masses (lbm)
    # batchVol:     batch volume at kegging (gal)
    #
    # OUTPUTS:
    # SRM:  predicted final color (SRM)
    
    # Pull lovibond color and mass info from grain bill
    L = grainBill['L']
    m = grainBill['m']
    
    # Calculate malt color units
    MCU = sum(L * m) / batchVol

    # Calculate SRM from MCU (morey)
    SRM = 1.4922 * MCU**0.6859
    
    # Output color prediction
    return(SRM)


def calcIBU(hopSchedule, SG_boil, boilVol):
    # Estimates batch IBUs using Tinseth's formulae from Palmer
    #
    # INPUTS:
    # hopSchedule:  dictionary containing
    #   m:  array of hop masses (oz)
    #   AA: array of alpha acid percentages (%)
    #   t:  array of hop addition times (min)
    # SG_boil:      pre-boil specific gravity
    # boilVol:      boil volume (gal)
    #
    # OUTPUTS:
    # IBU_batch:    Estimated batch bitterness (IBUs)
    
    # Pull mass, AA%, and boil time info from hop schedule
    m = hopSchedule['m']
    AA = hopSchedule['AA']
    t = hopSchedule['t']
    
    # Calculate array of Alpha Acid Units
    AAU = m * AA
    
    # Calculate time and SG utilization factors
    f_t = (1-np.exp(-0.04 * t)) / 4.15
    f_G = 1.65 * 0.000125**(SG_boil-1)

    # Calculate Utilization
    U = f_t * f_G

    # Calculate batch IBUs
    IBU_batch = np.sum(74.89 * AAU * U) / boilVol

    # Return IBUs
    return(IBU_batch)
    


# %% USE FUNCTIONS
# GOAL: build a gui for this as well

