# spect_lines
Spectral line finding algorithms and helper methods for the SNSPOL project. Made for the Hoffman group at the University of Denver

This is a work in progress project to create a spectral line finding algorithm for supernova spectroscopy although in theory it should work for any spectrum.

Pseudocode for the algorithm

maximums = empty list of values
minimums = empty list of values

for every flux value in spectrum_wavelengths:

  --if flux is higher than neighboring points:
    ----find left critical point
    ----find right critical point 
    ----if left & right critical points are far enough away from max:
      ------add wavelength and width values to maximums
  
  --if flux is lower than neighboring points:
    ----find left critical point
    ----find right critical point
    ----if left & right critical points are far enough away from min:
      ------add wavelength and width values to minimums
    
    
find_lines also contains the vel() method which calculates the relative velocities of detected line signatures around an inputted rest wavelength

There are many features to be added and changes to make, here's my TODO:

Add error handling
Implement error analysis and confidence metrics for line widths
Create new prominence checking method as an alternate filter method instead of line widths
Create epoch comparison method so that multiple epochs can be compared
Change critical point checking code so that slope is calculated once for the entire dataset
Reorganize code, relocating helper methods to be internal
General QOL changes and improving docstring wording
