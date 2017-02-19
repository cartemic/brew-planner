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
import re
import linecache


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


def infoGet(fName, key):
    # pulls a desired block of information from a recipe file
    #
    # INPUTS:
    #   fName:  string of the filepath for the desired recipe
    #   key:    string of the key for the desired information block, e.g.
    #           'grain' would return information about the grains etc.
    #
    # OUTPUT:
    # dictOut:  dictionary of desired information
    #
    # NOTE: this routine assumes that all attributes other than 'name' are
    #       floats and will convert them to a numpy float array without warning

    # open recipe file
    file = open(fName, 'r')

    # initialize blank dict for start/stop locations of desired info block
    locations = {}

    # search recipe for desired information block
    for num, line in enumerate(file, 1):
        if key+'-start' in line:
            locations['start'] = num+1
        elif key+'-end' in line:
            locations['end'] = num
            break

    # define proper start and end lines for block
    startLine = locations['start']
    endLine = locations['end']

    # close recipe file
    file.close()

    # initialize blank dictionary for desired output block
    dictOut = {}

    # pull information from recipe file
    for i in range(startLine, endLine):
        # pull lines and remove carriage return
        theLine = linecache.getline(fName, i)[:-1]

        # remove commas with regex
        pattern = re.compile(',')
        theLine = [x for x in pattern.split(theLine) if x]

        # build dictionary entry
        attribute = theLine[0]
        values = theLine[1:]
        if attribute != 'name':
            values = np.array(values, dtype='f')
        dictOut[attribute] = values

    # Return desired information
    return(dictOut)


def avgPitchRate(yeastType, OG):
    # Returns average yeast pitch rate depending on type of yeast and wort OG
    #
    # INPUTS:
    #   yeastType:  string of yeast type, either 'lager' or 'ale'
    #   OG:         wort OG
    #
    # OUTPUTS:
    #   pitchRate:  average pitch rate (million cells / ml - degP)
    #
    # NOTE:     any yeastType that isn't 'lager' or 'ale' will return NaN

    # estimates from NorthernBrewer yeast pitch pdf
    if yeastType == 'ale':
        if OG < 1.055:
            return(0.5)
        else:
            return(1.0)
    elif yeastType == 'lager':
        if OG < 1.055:
            return(1.0)
        else:
            return(1.5)
    else:
        return(np.NaN)


def yeastToPitch(pitchRate, OG, wortVol):
    # calculates the amount of yeast to pitch into the wort
    #
    # INPUTS:
    #   pitchRate:  desired pitch rate (million cells / ml - degP)
    #   OG:         wort OG
    #   wortVol:    wort volume (gal)
    #
    # OUTPUTS:
    #   pitchCell:  amount of yeast to pitch (billion cells)

    # convert to cells / ml - degP
    pitchRate = pitchRate * 1e6

    # convert OG to Plato
    OG = OG / 4

    # convert wort volume to ml
    wortVol = wortVol * 3785

    # calculate number of cells, convert to billions
    pitchCell = pitchRate * OG * wortVol / 1e9

    # return billions of cells
    return(pitchCell)


# %% USE FUNCTIONS
# GOAL: build a gui for this as well
