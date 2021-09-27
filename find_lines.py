# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 02:08:23 2021

@author: junec
"""

##Line signature location algorithm
##When inputting a spectrum will return a list of line signatures listed with global prominence 

import numpy as np
from scipy.constants import c


def find_lines(wave: list, flux: list, fit_length: int, min_width: float, wavelengths_only = False, weights : list = []):
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
        
        #print("Check width successfully called")
        
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
        
        Calculates the a least squares linear regression of the fluxes on both
        sides of a given index in order to find the critical points where the 
        line signature ends
    
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
    
    if wavelengths_only:
        maximums = [a[0] for a in maximums]
        minimums = [a[0] for a in minimums]
    
    return maximums, minimums
 

def vel(wave : list, flux : list, fit_length : int, min_width : float, rest_wavelength : float, wavelengths_only = False, weights=None):
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
    
    maximums, minimums = find_lines(wave, flux, fit_length, min_width, wavelengths_only, weights)
    
    def delta(wavelength : float):
        return ((c/1000) * (wavelength - rest_wavelength) / rest_wavelength)
    
    if wavelengths_only:
        maximums = [[a, delta(a)] for a in maximums if np.abs(delta(a)) < 30000]
        minimums = [[a, delta(a)] for a in minimums if np.abs(delta(a)) < 30000]
        return maximums, minimums
    
    maximums = [a + [delta(a[0])] for a in maximums if np.abs(delta(a[0])) < 30000]
    minimums = [a + [delta(a[0])] for a in minimums if np.abs(delta(a[0])) < 30000]
            
    return maximums, minimums

def validate_list(lists):
    
    if not isinstance(lists, list):
        return TypeError
    
    for n in lists:
        if not isinstance(n, (float, int)):
            return TypeError
        
    if len(lists) == 1:
        return ValueError
        
    return True