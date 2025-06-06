"""
The purpose of these sets of functions are to calculate the Zeta function of a time series, as described in 
Gomes et al. 
"""

#Imports:
from MFDFA import MFDFA
import matplotlib.pyplot as plt
import numpy as np

###############################################################################

# This is a function to calculate the lagged difference between B-field magnitude values
# (Essentially this performs the operation inside the square brackets in the above equation)
def lagged_difference(bmag, lag):
    """
    Calculates the lagged difference between B-field magnitude values
    
    Parameters
    ----------
    bmag : 1D Array
        Time series of B-field magnitude values

    lag : int
        The lag tau for a particular lagged difference
        
    Returns
    -------
    diff_list : 1D Array
        The list of differenced magnetic field magnitude         data
    """
    diff_list = []
    for t in range(len(bmag)-1):
        if t + lag <= len(bmag)-1:
            # Calculate the difference and append it to a              list:
            diff_list.append(np.abs(bmag[t+lag]-bmag[t]))
        else: break
    return np.array(diff_list)

###############################################################################

# This function will implement a calculation of the structure function:
def structure_function(bmag, plist,lags):
    """
    This function calculates the structure function of
    the magnetic field magnitude data. This represents the average magnitude of
    fluctuations of a signal at different scales. 

    Parameters
    ----------
    bmag : 1D Array
        Magnetic field magnitude time series.
    plist : List of Integers
        A list of scaling factors, p.
    lags : List of Integers
        List of integer lags. 

    Returns
    -------
    strucArray: 2D Array of structure functions for each scaling exponent p.
    """
    # The structure of this array will be [[S_p1(Tau_1), S_p1(Tau_2), S_p1(Tau_3),...],[S_p2(Tau_1), S_p2(Tau_2), S_p2(Tau_3),...],...]
    strucArray = np.zeros((len(plist),len(lags)))
    
    # Set up plot
    fig, ax = plt.subplots()

    # For each set of differenced lags, find the mean of 
    # the data to the power of p.
    for lag in lags:
        lag_diff = lagged_difference(bmag, lag)
        for p in range(len(plist)):
            strucArray[p-1][lag-1] = np.mean(lag_diff**p) 
    # Now plot the structure function for each p:
    for p in range(len(strucArray)):
        pArray = np.full(len(lags),plist[p])
        sc = ax.scatter(np.log(lags),np.log(strucArray[p]), c = pArray, s = 8, cmap = 'RdGy', vmin = -20,vmax = 20)    
    
    # Set the plot properties:
    cbar = plt.colorbar(sc)
    cbar.set_label('Order p')
    ax = plt.gca()
    ax.set_facecolor("black")
    plt.xlabel(r'$\log{\tau}$')
    plt.ylabel(r'$\log{S_p(\tau)}$')
    plt.title(r'$p$th-Order Structure Function $S_p(\tau)$ vs. Lag $\tau$ (Log-Scale)')
    plt.show()


    return strucArray

###############################################################################

def ZetaFunc(lags, strucArray):
    """
    Calculates the Zeta function, Zeta(p) from the input     structure function array. 

    Parameters
    ----------
    lags: 1D Array
        The set of lags used in the structure function.

    strucArray: 2D Array
        An array containing the structure function for           various lags and scaling exponents p.

    Returns
    -------
    zList : The set of Zeta(p) values calculated from the    structure function.
    """

    # Initialize an empty list to hold the Zeta function       values:
    zList = [] 
   
    # For every p in our scaling exponents, find the Zeta      Function values
    for struc in strucArray:
        # We fit to a 1st degree polynomial for a linear           fit
        Zeta = np.polyfit(np.log(lags), np.log(struc),1)[0]
        zList.append(Zeta)

    return zList 

###############################################################################

