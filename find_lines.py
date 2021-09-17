# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 02:08:23 2021

@author: junec
"""

##Line signature location algorithm
##When inputting a spectrum will return a list of line signatures listed with global prominence 

import numpy as np
from scipy.constants import c

dir = "C:/Users/junec/Documents/SN_2012au/"


def find_lines(wave: list, flux: list, fit_length: int, min_width: float, wavelengths_only = False):
    """
    Returns 2 two-dimensional lists listing the detected minima, maxima, and estimated line widths for the inputted spectrum
    
    Runs a line-finding algorithm which finds all local minima and maxima in the inputted spectrum
    Runs a second algorithm to estimate the line widths of the detected minima and maxima
    Filters out lines that are too narrow and creates 2 two-dimensional lists containing 
    the wavelength and index data for the detected minima and maxima

    Parameters
    ----------
    wave : list
        The wavelengths of the spectrum
    flux : list
        The fluxes corresponding to the inputted wavelengths
    fit_length : int
        The length of the linear regression in indices used to find the estimated width of each line
    min_width : float
        The minimum acceptable width in angstroms for detected lines
    wavelengths_only : boolean
        If set to 'True' returns only the wavelengths of detected line signatures

    Returns
    -------
    maximums : list
        A list containing the wavelength and line widths of each detected maximum with each entry formatted as:
        [wavelength, left_width_Å, right_width_Å, wavelength_indices, left_width_indices, right_width_indices]
    minimums : list
        A list containing the wavelength and line widths of each detected minimum with each entry formatted as:
        [wavelength, left_width_Å, right_width_Å, wavelength_indices, left_width_indices, right_width_indices]

    """
    
    maximums = [[0,0,0,0,0,0]]
    minimums = [[0,0,0,0,0,0]]
    
    λ = round(wave[1] - wave[0])
    min_width = min_width / λ
    
    # checks each wavelength for maxima and minima
    
    for index in range(len(wave)-1):
          
          # checking for maximum
          # maximum needs to have a higher flux than the neighboring wavelengths
          
          if flux[index] > flux[index-1] and flux[index] > flux[index+1]:
              
              left, right = check_width(wave, flux, index, fit_length, 1)
              
              # filters out narrow lines
              
              if(left > min_width and right > min_width):
                  
                  if(maximums == [[0,0,0,0,0,0]]): # special case for initializing the list
                  
                      maximums = [[wave[index], left*λ, right*λ, index, left, right]]
                  
                  else: 
                      
                      maximums.append([wave[index], left*λ, right*λ, index, left, right])
           
          # checking for minimum
          # minimum needs to have a lower flux than the neighboring wavelengths  
              
          if flux[index] < flux[index-1] and flux[index] < flux[index+1]:
              
              left, right = check_width(wave, flux, index, fit_length, -1)
              
              # filters out narrow lines
              
              if(left > min_width and right > min_width):
                  
                   if(minimums == [[0,0,0,0,0,0]]): # special case for initializing the list
                   
                      minimums = [[wave[index], left*λ, right*λ, index, left, right]]
                  
                   minimums.append([wave[index], left*λ, right*λ, index, left, right])
                   
    if wavelengths_only:
        maximums = [a[0] for a in maximums]
        minimums = [a[0] for a in minimums]
    
    return maximums, minimums
    

def check_width(wave : list, flux : list, i : int, fit_length : int, maxmin: int):
    """
    Internal helper method which returns the width of an inputted line
    
    Checks the width of an inputted line signature by fitting a line to the flux
    signature on either side of the minimum or maximum and checking if fitted
    line has hit a critical point. The estimated location of that critical point
    is returned as the width of the line signature on that side

    Parameters
    ----------
    wave : list
        The wavelengths of the spectrum.
    flux : list
        The fluxes corresponding to the inputted wavelengths
    i : int
        The index of the detected minimum or maximum
    fit_length : int
        The length of the linear regression in indices used to find the estimated width of each line
    maxmin : int
        Changes the mode of check_width. maxmin = -1 yields minimum mode, maxmin = +1 yields maximum mode

    Returns
    -------
    left_width : float
        The left width of the inputted line signature
    right_width : float
        The right width of the inputted line signature

    """
    
    #print("Check width successfully called")
    
    index = i
    left_width = 0
    right_width = 0
    
    # maximum mode
    
    if maxmin == 1:       
        
        # checking for left critical point 
        
        while index != -1:
            
            index = index - 1
            
# =============================================================================
#             # Line width algorithm that used averages to find line widths
#             # Currently unused
#
#             left_flux, right_flux = check_average(wave, flux, i, fit_length)
#             
#             if left_flux > flux[index+1]:
#                 left_width = i - index
#                 break
#             #if no point to the left of the peak is higher than peak then left_width found
#             if index == 0:
#                 left_width = i
# =============================================================================

            left_slope, right_slope = ls_line(wave, flux, index, fit_length)
            
            if left_slope < 0:
                left_width = i - (index + fit_length/2)
                break
            if index == 0: # catches special case where there are no wavelengths on left side
                left_width = i
        
        index = i
        
        # checking for right critical point
        
        while index != len(flux):
            
            index = index + 1 
            
# =============================================================================
#             # Line width algorithm that used averages to find line widths
#             # Currently unused
#
#             left_flux, right_flux = check_average(wave, flux, i, fit_length)
#             
#             if right_flux > flux[index-1]:
#                 right_width = index - i
#                 break
#             #if no point to the right of the peak is higher than peak then right_width found
#             if index == len(flux) - 1:
#                 right_width = i
# =============================================================================

            left_slope, right_slope = ls_line(wave, flux, index, fit_length)
            
            if right_slope > 0:
                right_width = (index + fit_length/2) - i
                break
            if index == len(flux) - 1: # catches special case where there are no wavelengths on right side
                right_width = i
               
    # minimum mode
    
    if maxmin == -1:

        # checking for left critical point         

        while index != -1:
            
            index = index - 1
            
# =============================================================================
#             # Line width algorithm that used averages to find line widths
#             # Currently unused
#
#             left_flux, right_flux  = check_average(wave, flux, i, fit_length)
#             
#             if left_flux < flux[index+1]:
#                 left_width = i - index
#                 break
#             if index == 0:
#                 left_width = i
# =============================================================================
            
            left_slope, right_slope = ls_line(wave, flux, index, fit_length)
            
            if left_slope > 0:
                left_width = i - (index + fit_length/2)
                break
            if index == 0: # catches special case where there are no wavelengths on left side
                left_width = i
        
        index = i
        
        # checking for right critical point
        
        while index != len(flux):
            
            index = index + 1
            
# =============================================================================
#             # Line width algorithm that used averages to find line widths
#             # Currently unused
#
#             left_flux, right_flux = check_average(wave, flux, i, fit_length)
#             
#             if right_flux < flux[index-1]:
#                 right_width = index - i
#                 break
#             if index == len(flux) - 1:
#                 right_width = i
# =============================================================================
            
            left_slope, right_slope = ls_line(wave, flux, index, fit_length)
            
            if right_slope < 0:
                right_width = (index + fit_length/2) - i
                break
            if index == len(flux) - 1: # catches special case where there are no wavelengths on right side
                right_width = i
                
                
    return left_width, right_width


def check_average(wave : list, flux : list, i : int, fit_length : int):
    """
    Returns the arithmetic mean of the fluxes on either side of an inputted wavelength
    
    ***CURRENTLY UNUSED***

    Parameters
    ----------
    wave : list
        The wavelengths of the spectrum.
    flux : list
        The fluxes corresponding to the inputted wavelengths
    i : int
        The index of the detected minimum or maximum
    fit_length : int
        The length of the arithmetic mean in indices used to find the estimated width of each line

    Returns
    -------
    left_average : float
        The arithmetic mean of the fluxes on the left side of i
    right_average : float
        The arithmetic mean of the fluxes on the right side of i

    """
    
    left_average = 0
    right_average = 0
    
    # uses counters to deal with edge case where investigated average is near to the edge of the dataset
    
    left_length = 0
    right_length = 0
    
    for n in range(fit_length):
        if i - (n+1) > 0:
            left_average = left_average + flux[i - (n+1)]
            left_length = left_length + 1
        if i + (n+1) < len(flux):    
            right_average = right_average + flux[i + (n+1)]
            right_length = right_length + 1
    
    if left_length > 0 and right_length > 0:
        left_average = left_average / left_length
        right_average = right_average / right_length
    
    return left_average, right_average


def ls_line(wave : list, flux : list, i : int, fit_length : int):
    """
    Helper function for check_width, returns the slope of the data on either side of an inputted wavelength
    
    Calculates the a least squares linear regression of the fluxes on both
    sides of a given index in order to find the critical points where the 
    line signature ends

    Parameters
    ----------
    wave : list
        The wavelengths of the spectrum.
    flux : list
        The fluxes corresponding to the inputted wavelengths
    i : int
        The index of the detected minimum or maximum
     fit_length : int
        The length of the linear regression in indices used to find the estimated width of each line

    Returns
    -------
    left_slope : flaot
        The calculated slope to the left of the given wavelength
    right_slope : float
        The calculated slope to the right of the given wavelength

    """
    
# =============================================================================
#     Least square linear regression assumes that there is a linear mathematical relationships between two discrete variables
#     This routine calculates the slope of this theoretical line as a way of finding where the slope changes signs
#     Slope changing signs indicates the end of a maximum or minimum and is used as the basis for determining line width
#
#     Assumes data takes form y = A + Bx
#     A = (Σx^2 * Σy - Σx * Σxy) / Δ
#     B = (N * Σxy - Σx * Σy) / Δ
#     Δ = N * Σx^2 - (Σx)^2
#
#     N is the number of data points
#
#     As this routine only needs B to find slope inflection A is not calculated
#     Slope is calculated on each side of given index out to the given width
# =============================================================================
    
    left_xy = 0
    left_y = 0
    left_x = 0
    left_x2 = 0
    
    right_xy = 0
    right_y = 0
    right_x = 0
    right_x2 = 0
    
    for n in range(fit_length):
        if i - (n+1) > 0:
            left_xy = left_xy + flux[i - (n+1)] * wave[i - (n+1)]
            left_y = left_y + flux[i - (n+1)]
            left_x = left_x + wave[i - (n+1)]
            left_x2 = left_x2 + np.power(wave[i - (n+1)],2)
        if i + (n+1) < len(flux):
            right_xy = right_xy + flux[i + (n+1)] * wave[i + (n+1)]
            right_y = right_y + flux[i + (n+1)]
            right_x = right_x + wave[i + (n+1)]
            right_x2 = right_x2 + np.power(wave[i + (n+1)],2)
    
    left_slope = (fit_length * left_xy - left_x * left_y) / (fit_length * left_x2 - np.power(left_x, 2))
    right_slope = (fit_length * right_xy - right_x * right_y) / (fit_length * right_x2 - np.power(right_x, 2))
    
    return left_slope, right_slope
   

def vel(wave : list, flux : list, fit_length : int, min_width : float, rest_wavelength : float, max_velocity : float):
    """
    Finds the relative velocities of all line signatures around an inputted rest wavelength

    Parameters
    ----------
    wave : list
        The wavelengths of the spectrum
    flux : list
        The fluxes corresponding to the inputted wavelengths
    fit_length : int
        The length of the linear regression in indices used to find the estimated width of each line
    min_width : float
        The minimum acceptable width in angstroms for detected lines
    rest_wavelength : float
        The wavelength around which all the relative velocities are calculated
    max_velocity : float
        The maximum relative velocity allowed

    Returns
    -------
    maximums : list
        A list containing the wavelength, line widths, and relative velocities
        of each detected maximum with each entry formatted as:
        [wavelength, left_width_Å, right_width_Å, wavelength_indices, left_width_indices, right_width_indices, relative_velocity]
    minimums : list
        A list containing the wavelength, line widths, and relative velocities
        of each detected minimum with each entry formatted as:
        [wavelength, left_width_Å, right_width_Å, wavelength_indices, left_width_indices, right_width_indices, relative_velocity]

    """
    
    maximums, minimums = find_lines(wave, flux, fit_length, min_width)

    maximums = [a + [((c/1000) * (a[0] - rest_wavelength) / rest_wavelength)] for a in maximums if ((c/1000) * np.abs(a[0] - rest_wavelength) / rest_wavelength) < max_velocity]
    minimums = [a + [((c/1000) * (a[0] - rest_wavelength) / rest_wavelength)] for a in minimums if ((c/1000) * np.abs(a[0] - rest_wavelength) / rest_wavelength) < max_velocity]
            
    return maximums, minimums