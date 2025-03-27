"""
  This is a library that utilizes the MFDFA and Pyspedas libraries in order to
  perform Multifractal Detrended Fluctuation Analysis on solar wind magnetic
  field data in sub and super Alfvenic regions 

Created by : Sam Coscia

"""
#Imports:
from MFDFA import MFDFA
import pyspedas
import matplotlib.pyplot as plt
import numpy as np
from pyspedas import tplot
from scipy import signal 

#############################################################################

def MFDFA_analysis(timeSeries, timeSeries2 = None):

    """
    This is a function that performs the entire MFDFA process on an input time
    series. It can also take in a second time series, and plot their spectra against
    each other

    Parameters
    ----------
    timeSeries : 1D array
        The first magnetic field time series data to perform MFDFA on.
    timeSeries2 : 1D array
        An optional second magnetic field time series data entry to simultaneously
        perform MFDFA on.
    Returns
    ------ 
    None (Only plots)

    """

    # We first need to select a band of lags, varying from small segments of data to long ones:
    # These bands must be integers, since they will segment the data into chunks of side length s:
    lag = np.arange(1,len(timeSeries),1000)

    # We need to select a range of powers for the fluctuation function.
    # Gomes et al. do -20 to 20 with increments of 0.25
    qList=np.arange(-20,20,0.25).tolist()
    qList.remove(0)
    
    #The order for the polynomial fitting is chosen by Gomes et al. as 3:
    order=3
    
    # If only 1 time series is given, only produce plots for 1 time series : 
    if timeSeries2 is None:
    # Find the fluctuation function for each q, and save them to a list:
        lag,dfaList=flucFunc(timeSeries,lag,qList,order)

    # Here I am extracting the Hurst Parameter for each fluctuation function.
    # I use the transpose of dfaList for the purpose of easily parsing each row
    # of the matrix in the form : [[F_q1(s1), F_q1(s2), ...], [F_q2(s1), F_q2(s2), ...], ...]
        h=hList(lag, qList, dfaList.T)
        ########PLOTTING#################
        plt.plot(qList,h)
        plt.ylabel('h(q)')
        plt.xlabel('q')
        plt.title('Hurst Exponent Spectrum')
        plt.show()
        #################################

        # Now find and plot the Renyi exponent spectrum:
        tauFunc=renyiExp(qList,h)
        ########PLOTTING#################
        plt.plot(qList,tauFunc,'o')
        plt.xlabel(r'$q$')
        plt.ylabel(r'$\tau(q)$')
        plt.title('Renyi Exponent Spectrum')
        plt.show()
        #################################

        # Finally, from the Renyi exponents, find and plot the multifractal spectrum:
        alphaList, mfList = mfSpec(tauFunc, qList)
        #######PLOTTING##################
        # Plot the MF spectrum function f(alpha)
        plt.plot(alphaList,mfList)
        plt.title('Singularity Spectrum')
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$f(\alpha)$')
        plt.show()
        #################################
        # Here we are calculating the PSD function:
        # Choose 22 Hz as the sampling frequency, as in Gomes et al.
        fs = 22 # Hz
        fList, PSDlist = PSDfunc(timeSeries,fs) 
        #######PLOTTING###################
        plt.loglog(f,PSD)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Power Spectral Density [nt**2/Hz]')  
        ##################################

    # If 2 time series are given, perform the MFDFA process on both, and
    # superimpose the plots:
    else:
        
        # Fluctuation functions :  
        lag1,dfaList1 = flucFunc(timeSeries,lag,qList,order)
        lag2,dfaList2 = flucFunc(timeSeries2,lag,qList,order)
        
		# Hurst parameters from the fluctuation function:
        h1=hList(lag, qList, dfaList1.T)
        h2=hList(lag, qList, dfaList2.T)
        ########PLOTTING#################
        plt.plot(qList,h1,label = "Sub-Alfvenic Interval (April 28th, 2021 from 09:33-14:42 UTC)")
        plt.plot(qList,h2,label = "Super-Alfvenic Interval (April 28th, 2021 from 00:00-09:33 UTC")
        plt.ylabel('h(q)')
        plt.xlabel('q')
        plt.title('Hurst Exponent Spectrum')
        plt.show()
        #################################

        # Now find and plot the Renyi exponent spectrum:
        tauFunc1=renyiExp(qList,h1)
        tauFunc2=renyiExp(qList,h2)
        ########PLOTTING#################
        plt.plot(qList,tauFunc1,'o',label = 'Sub-Alfvenic Interval (April 28th, 2021 from 09:33-14:42 UTC)')
        plt.plot(qList,tauFunc2,'o', label = 'Super-Alfvenic Interval (April 28th, 2021 from 00:00-09:33 UTC)')
        plt.xlabel(r'$q$')
        plt.ylabel(r'$\tau(q)$')
        plt.title('Renyi Exponent Spectrum')
        plt.show()
        #################################

        # Finally, from the Renyi exponents, find and plot the multifractal spectrum:
        alphaList1, mfList1 = mfSpec(tauFunc1, qList)
        alphaList2, mfList2 = mfSpec(tauFunc2, qList)
        #######PLOTTING##################
        # Finally, plot the f(\alpha) function
        plt.plot(alphaList1, mfList1, label = 'Sub-Alfvenic Interval (April 28th, 2021 from 09:33-14:42 UTC)')
        plt.plot(alphaList2, mfList2, label = 'Super-Alfvenic Interval (April 28th, 2021 from 00:00-09:33 UTC)')
        plt.title('Singularity Spectrum')
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$f(\alpha)$')
        plt.show()
        #################################

        # Here we are calculating the PSD function:
        # Choose 22 Hz as the sampling frequency, as in Gomes et al.
        fs = 22 # Hz
        fList1, PSDlist1 = PSDfunc(timeSeries1,fs) 
        fList2, PSDlist2 = PSDfunc(timeSeries2,fs) 
        #######PLOTTING###################
        # Here we plot the PSD function on a log-log scaling:
        plt.loglog(fList1,PSDlist1, label = 'Sub-Alfvenic Interval (April 28th, 2021 from 09:33-14:42 UTC)')
        plt.loglog(fList2,PSDlist2, label = 'Super-Alfvenic Interval (April 28th, 2021 from 00:00-09:33 UTC)')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Power Spectral Density [nt**2/Hz]')  
        ##################################

    return None

#############################################################################

def flucFunc(bmag, lag, qList, order):
    """  
    In this function, we calculate the fluctuation function of input magnetic
    field time series data.    
    
    Parameters
    ---------- 
    bmag : 1D array
       The magnetic field time series.
    lag : 1D array
       An array containing increasing integer values for each lag, s in the 
       fluctuation function.
    qList : 1D array
       An array containing the scaling powers q to evaluate the fluctuation
       function at.
    order : int
       The order of the polynomial fitting performed in the calculation of the
       fluctuation function. 

    Returns 
    -------
    lag : 1D array 
       The list of lags used in the fluctuation function.
    dfaList : 2D array
       The fluctuation function evaluated at each value of q and for each lag.
         """
    # Here we use the MFDFA function to give us the fluctuation function ,
    # which we then take the log of:
    lag, dfaList = MFDFA(bmag, lag = lag, q = qList, order = order)

    # Initialize the plots: 
    fig, ax = plt.subplots()

    """
       NOTE : 
       The structure of dfaList is of the form: [[F_s1(q1), F_s1(q2),...],[F_s2(q1),
       F_s2(q2), ...], ...], but to plot the fluctuation function using a different 
       color for each value of q, we take the transpose of dfaList, to get an array 
       of the form : 
           [[F_q1(s1), F_q1(s2), ...], [F_q2(s1), F_q2(s2), ...], ...].
     
    """
    # First determine the log of the lags and the fluctuation matrix
    logdfa = np.log(dfaList.T)
    logs = np.log(lag)

    # for every q, determine the F_q vs. s plot.
    for q in range(len(dfaList.T)): # The transpose is taken here!
        # Here we plot the fluctuation function with a color bar for each q value:
        qArray=np.full(len(lag),qList[q])
        sc = ax.scatter(logs,logdfa[q],c = qArray,s = 8, cmap ='inferno', vmin=-20, vmax=20)
    
    # This is just some plot style related code:
    cbar = plt.colorbar(sc)
    cbar.set_label('Order q')
    ax = plt.gca()
    ax.set_facecolor("black")
    plt.xlabel(r'$\log{s}$')
    plt.ylabel(r'$\log{F_q(s)}$')
    plt.title(r'$q$th-Order Fluctuation Function $F_q(s)$ vs. Lag $s$ (Log-Scale)')
    plt.show()
    
    return lag,dfaList

#############################################################################

def hList(lag, dfaListT):
    """
    Calculates the Hurst spectrum, h(q) from the input fluctuation function
    data. 

    Parameters
    ----------
    lag : 1D array
        The set of lags used in the fluctuation function.
    dfaListT : 2D array
        The transpose of the fluctuation function evaluated for each lag and
each value of q in the form : [[F_q1(s1),F_q1(s2),...],[F_q2(s1),F_q2(s1),...]]
    qList : 1D array 
        An array containing the scaling exponents of the fluctuation function,
        q. 
    Returns
    -------         
    hList : 1D array
        The Hurst spectra, h(q), over all q values.
 
    """
    # Initialize an empty list to hold our Hurst parameters:
    hList=[]
    
    # For every q in our list of powers, we find the Hurst parameter:
    for dfa in dfaListT:
        # We fit to a 1st degree polynomial for a linear fit
        H = np.polyfit(np.log(lag[2:len(lag)]),np.log(dfa[2:len(lag)]),1)[0]
        hList.append(H)

    return hList

#############################################################################

def renyiExp(qList,hList):
    """
    This function calculates the Renyi exponents from the Hurst spectrum, using
    the relation: tau(q)=q*h(q)-1. 
    
    Parameters
    ----------
    qList : 1D array 
        An array containing the scaling exponents of the fluctuation function,
        q. 
    hList : 1D array
        An array containing the hurst spectrum values for each q.

    Returns
    -------
    tauList : 1D array
        An array containing the Renyi exponent spectrum values for each q. 

    """

    #Initialize an empty list to hold the Renyi exponents:
    tauList=[]
    
    # Now, for each q and H, we can calculate the Renyi exponent defined by
    # the above equation 
    for i in range(len(qList)):
        # Here we find the Renyi exponent and append it to a list :
        tauList.append(qList[i]*hList[i]-1)
    
    return tauList        

#############################################################################

def mfSpec(tauFunc,qList):
    """
    In this function we calculate the multifractal spectrum from the Renyi
	exponent spectrum using the equation:

                 f(alpha) = q*alpha-tau(q), where alpha = tau'(q)
   
	Parameters     
	----------
    tauFunc : 1D array
        An array containing the Renyi exponent values for each q.
    qList : 1D array 
        An array containing the scaling exponents of the fluctuation function,
        q. 
 
    Returns
	-------
    alphaList : 1D array
        An array containing each of the values of alpha = tau'(q) for every q. 

    mfList : 1D array
        An array containing each of the values of the MF spectrum, f(alpha)

    """    
    # First we take the derivative of tau(q) to find alpha:
    alphaList=np.gradient(tauFunc,qList)

    #Initialize an empty list that will hold the f(alpha) values:
    mfList=[]
    for i in range(len(tauFunc)):
        #Perform the operation for f(alpha) written above:
        mfList.append(qList[i]*alphaList[i]-tauFunc[i])

    return alphaList, mfList

#############################################################################

def PSDfunc(bmag, fs):
    """
    This is a function to perform all Power Spectral Density related
    calculations.

    Parameters
    ----------
    bmag : 1D array
       An array containing the magnetic field time series data.
    fs : int
       The sampling frequency for calculating the power spectral density.

    Returns
    ------- 
    f : 1D array
       An array of frequencies the PSD function is evaluated at.
    PSD : 1D array
       The values of the PSD function for each frequency. 
    """

    # The signal.periodogram function calculates the Power Spectral Density for
    # us :
    f, PSD = signal.periodogram(bmag,fs)
  
    return f,PSD

#############################################################################

def shuffle(timeSeries):
    """
    This function randomly shuffles an input time series.

    Parameters
	----------
	timeSeries : 1D array
       The input time series data.   	
    
	Returns
	-------
    to_shuffle : 1D array
       The shuffled time series data. 
    """

    # First we copy the input time series:
    to_shuffle = np.copy(timeSeries)
    # Use the default rng generator:
    rng = np.random.default_rng()
    # Finally, shuffle the copied time series using the selected rng generator:
    rng.shuffle(to_shuffle)

    return to_shuffle

#############################################################################

def magnitude(list):
    """
    This function simply calculates the magnitude of vectors in a list. 
    PSP magetic field data comes only in RTN (radial-tangential-normal)
	coordinates. We can use this function to calculates the magnitude of
    the magnetic field from the RTN coordinates. 

    Parameters
    ----------
    list : list of vectors
       The input vector list.

    Returns
	-------
    bmagList : 1D array 
       The magnetic field magnitude data.

    """
    bmagList = []

    for vector in list:
        # Calculate the magnitude here:
        bmag = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)  
        bmagList.append(bmag)
  
    return np.array(bmagList)

#############################################################################

def phaseShuffle(timeSeries):
    """
    This is a function to randomize the phases of a time series.
    
    Parameters
    ----------
    timeSeries : 1D array
        Input time series data f(t).

    Returns
    -------
    shuffled_timeSeries : 1D array
        The input time series with randomly shuffled phases.
    """    

    # First we take the FFT of the time series
    FFT = np.fft.fft(timeSeries)
    
    # Now find the phases and magnitudes : 
    phases = np.angle(FFT)
    mag = np.abs(FFT)

    # Randomize the phases: 
    np.random.shuffle(phases)

    # Now reconstruct the Fourier transform with randomized phases
    FFT_shuffled = mag * np.exp(1j * phases)

    # Finally, reconstruct the original function with the IFT
    shuffled_timeSeries = np.fft.ifft(FFT_shuffled)

    return shuffled_timeSeries
