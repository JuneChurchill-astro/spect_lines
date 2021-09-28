# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 02:08:23 2021

@author: junec
"""

## Line signature location algorithm
## When inputting a spectrum will return a list of line signatures listed with estimated widths

import numpy as np
from scipy.constants import c


def find_lines(wave: list, flux: list, fit_length: int, min_width: float, wavelengths_only = False, weights : list = [], com = True):
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
        If True returns maximums and minimums with only wavelengths listed, default : False
    weights : list
        The weights for the weighted linear regression
    con : boolean
        If true consolidates overlapping lines

    Returns
    -------
    maximums : list
        A list containing the wavelength and line widths of each detected maximum with each entry formatted as:
        [wavelength, left_width_Å, right_width_Å, wavelength_indices, left_width_indices, right_width_indices]
    minimums : list
        A list containing the wavelength and line widths of each detected minimum with each entry formatted as:
        [wavelength, left_width_Å, right_width_Å, wavelength_indices, left_width_indices, right_width_indices]

    """
    
    def check_width(i : int, maxmin: int):
        """
        Internal helper method which returns the width of an inputted line
        
        Checks the width of an inputted line signature by fitting a line to the flux
        signature on either side of the minimum or maximum and checking if fitted
        line has hit a critical point. The estimated location of that critical point
        is returned as the width of the line signature on that side
    
        Parameters
        ----------
        i : int
            The index of the detected minimum or maximum
        maxmin : int
            Changes the mode of check_width. maxmin = -1 yields minimum mode, maxmin = +1 yields maximum mode
    
        Returns
        -------
        left_width : float
            The left width of the inputted line signature
        right_width : float
            The right width of the inputted line signature
    
        """
        
        index = i
        left_width = 0
        right_width = 0
        
        # maximum mode
        
        if maxmin == 1:       
            
            # checking for left critical point 
            
            while index != -1:
                
                index = index - 1

                left_slope = ls_line(index, -1)
                
                if left_slope < 0:
                    
                    left_width = i - (index + fit_length/2)
                    break
                
                if index == 0: # catches special case where there are no wavelengths on left side
                
                    left_width = i
            
            index = i
            
            # checking for right critical point
            
            while index != len(flux):
                
                index = index + 1 
    
                right_slope = ls_line(index, 1)
                
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
                
                left_slope = ls_line(index, -1)
                
                if left_slope > 0:
                    left_width = i - (index + fit_length/2)
                    break
                if index == 0: # catches special case where there are no wavelengths on left side
                    left_width = i
            
            index = i
            
            # checking for right critical point
            
            while index != len(flux):
                
                index = index + 1
                
                right_slope = ls_line(index, 1)
                
                if right_slope < 0:
                    
                    right_width = (index + fit_length/2) - i
                    break
                
                if index == len(flux) - 1: # catches special case where there are no wavelengths on right side
                
                    right_width = i
                    
                    
        return left_width, right_width
    
    
    def ls_line(i : int, direction : int):
        """
        Helper function for check_width, returns the slope of the data on either side of an inputted wavelength
        
        Calculates a weighted least squares linear regression of the fluxes on both
        sides of a given index. Finds the critical points where the line signature ends
    
        Parameters
        ----------
        i : int
            The index of the detected minimum or maximum
        direction : int
            Determines which side the fit is calculated for, -1 = 'left' and 1 = 'right'
    
        Returns
        -------
        w_wlr[1] : float
            The slope of the calculated linear regression
    
        """
        
        X0, y, res = np.zeros(fit_length), np.zeros(fit_length), np.ones(fit_length)
        
        for n in range(fit_length):
            
            if i - (n+1) > 0 and direction == -1:
                
                X0[n] = wave[i - (fit_length - n)]
                y[n] = flux[i - (fit_length - n)]
                
                if np.any(weights) and len(weights) == len(wave):
                    
                   res[n] = weights[i - (fit_length - n)]
            
            if i + (n+1) < len(flux) and direction == 1:
                
                X0[n] = wave[i + n]
                y[n] = flux[i + n]
                
                if np.any(weights) and len(weights) == len(wave):
                    
                    res[n] = weights[i + n]
                    
        C = np.diag(res**2)
        
        #MLE least squares weighted linear regression
        
        X = np.c_[np.ones(X0.shape),X0]
        
        w_wlr = np.linalg.pinv(X.T @ np.linalg.pinv(C) @ X) @ (X.T @ np.linalg.pinv(C) @ y)
        
        return w_wlr[1]
    
    
    def combine(lines : list, maxmin : int, i = 0):
        """
        Combines overlapping minima and maxima
        
        Parameters
        ----------
        lines : list
            The minima or maxima that need to be checked for overlapping signatures
        maxmin : int
            Determines the mode, 1 = max, -1 = min
        
        Returns
        -------
        lines : list
            The combined signatures
        """

        if len(lines) > 1 and i < len(lines) - 1:
            
                ## if the bounds of one line overlap the center of the other and vice versa they are overlapping
                ## combines overlapping lines into one signature
                ## sets new center line of combined signature to be the local extrema of the two combining lines
                ## cleans up maxima and minima from noisy data
                
                if (lines[i][0] + lines[i][2] >= lines[i+1][0] and lines[i+1][0] - lines[i+1][1] <= lines[i][0]): 
                    
                    # checking minima
                    
                    if maxmin == -1:
                        
                        lines[i+1] = (lines[i][0], lines[i][1], lines[i+1][0]-lines[i][0]+lines[i+1][2],
                                      lines[i][3], lines[i][4], lines[i+1][3]-lines[i][3]+lines[i+1][5]
                                      ) if flux[lines[i][3]] <= flux[lines[i+1][3]] else ( # finds local minimum of two lines
                                          lines[i+1][0], lines[i+1][0]-lines[i][0]+lines[i][1], lines[i+1][2],
                                          lines[i+1][3], lines[i+1][3]-lines[i][3]+lines[i][4], lines[i+1][5])
                    
                    # checking maxima
                                          
                    if maxmin == 1:    
                        
                        lines[i+1] = (lines[i][0], lines[i][1], lines[i+1][0]-lines[i][0]+lines[i+1][2],
                                      lines[i][3], lines[i][4], lines[i+1][3]-lines[i][3]+lines[i+1][5]
                                      ) if flux[lines[i][3]] >= flux[lines[i+1][3]] else ( # finds local maximum of two lines
                                          lines[i+1][0], lines[i+1][0]-lines[i][0]+lines[i][1], lines[i+1][2],
                                          lines[i+1][3], lines[i+1][3]-lines[i][3]+lines[i][4], lines[i+1][5])
                    
                if lines[i+1][0] == lines[i][0]: # checks if two lines were just combined
                
                    del(lines[i]) # deletes duplicate line
                    lines = combine(lines, maxmin, i) # checks if any other lines lie within the new combined line
                    
                else: 
                    lines = combine(lines, maxmin, i + 1) # if no other lines overlapping on current index moves to next line
        
        return lines

    
    maximums = [[0,0,0,0,0,0]]
    minimums = [[0,0,0,0,0,0]]
    
    λ = round(wave[1] - wave[0])
    min_width = min_width / λ
    
    try:
        validate_list(wave)
    except TypeError:
        raise TypeError('Wavelengths not inputted as a list of floats or ints')
    except ValueError:
        raise ValueError('Inputted lists must be at least of length 2')
    
    # checks each wavelength for maxima and minima
    
    for index in range(len(wave)-1):
          
          # checking for maximum
          # maximum needs to have a higher flux than the neighboring wavelengths
          
          if flux[index] > flux[index-1] and flux[index] > flux[index+1]:
              
              left, right = check_width(index, 1)
              
              # filters out narrow lines
              
              if(left > min_width and right > min_width):
                  
                  if(maximums == [[0,0,0,0,0,0]]): # special case for initializing the list
                  
                      maximums = [[wave[index], left*λ, right*λ, index, left, right]]
                  
                  else:
                      
                      maximums.append([wave[index], left*λ, right*λ, index, left, right])
           
          # checking for minimum
          # minimum needs to have a lower flux than the neighboring wavelengths  
              
          if flux[index] < flux[index-1] and flux[index] < flux[index+1]:
              
              left, right = check_width(index, -1)
              
              # filters out narrow lines
              
              if(left > min_width and right > min_width):
                  
                   if(minimums == [[0,0,0,0,0,0]]): # special case for initializing the list
                   
                      minimums = [[wave[index], left*λ, right*λ, index, left, right]]
                  
                   else: 
                       
                      minimums.append([wave[index], left*λ, right*λ, index, left, right])
    
    # combines overlapping lines
    
    if com:
        
        maximums = combine(maximums, 1)
        minimums = combine(minimums, -1)
    
    if wavelengths_only:
        
        maximums = [a[0] for a in maximums]
        minimums = [a[0] for a in minimums]
    
    return maximums, minimums
 

def vel(wave : list, flux : list, fit_length : int, min_width : float, rest_wavelength : float, wavelengths_only = False, weights = None, com = True):
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
    wavelengths_only : boolean
        If True returns maximums and minimums with only wavelengths listed, default : False
    weights : list
        The weights for the weighted linear regression

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
    
    maximums, minimums = find_lines(wave, flux, fit_length, min_width, wavelengths_only, weights, com)
    
    def delta(wavelength : float): # calculates line velocity from rest_wavelength
        return ((c/1000) * (wavelength - rest_wavelength) / rest_wavelength)
    
    if wavelengths_only:
        
        maximums = [[a, delta(a)] for a in maximums if np.abs(delta(a)) < 30000]
        minimums = [[a, delta(a)] for a in minimums if np.abs(delta(a)) < 30000]
        
        return maximums, minimums
    
    maximums = [a + [delta(a[0])] for a in maximums if np.abs(delta(a[0])) < 30000]
    minimums = [a + [delta(a[0])] for a in minimums if np.abs(delta(a[0])) < 30000]
            
    return maximums, minimums


def validate_list(lists):
    """
    Checks if user input is of a valid format to be analyzed
    
    Parameters
    ----------
    lists : list
        User inputted data either wavelengths or flux
        
    Returns
    -------
    TypeError
        If data is not a list or is empty list
    ValueError
        If data is only 1 entry long
    True
        If data is valid for analysis
    """
    
    if not isinstance(lists, list):
        
        return TypeError
    
    for n in lists:
        
        if not isinstance(n, (float, int)):
            
            return TypeError
        
    if len(lists) == 1:
        
        return ValueError
        
    return True
