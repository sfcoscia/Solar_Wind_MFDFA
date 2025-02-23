# This is a library that utilizes the MFDFA and Pyspedas libraries in order to
# perform Multifractal Detrended Fluctuation Analysis on solar wind magnetic
# field data in sub and super Alfvenic regions 

#Imports:
from MFDFA import MFDFA
import pyspedas
import matplotlib.pyplot as plt
import numpy as np
from pyspedas import tplot
from scipy import signal 

#############################################################################

#  This is the main function we will use to run our helper functions:

def MFDFA_analysis(timeSeries):
    # We first need to select a band of lags, varying from small segments of data to long ones:
    # These bands must be integers, since they will segment the data into chunks of side length s:
    lag = np.arange(1,len(timeSeries),1000)

    # We need to select a range of powers for the fluctuation function.
    # Gomes et al. do -20 to 20 with increments of 0.25
    qList=np.arange(-20,20,0.25).tolist()
    qList.remove(0)
    
    #The order for the polynomial fitting is chosen by Gomes et al. as 3:
    order=3
    # Find the fluctuation function for each q, and save them to a list:
    lag,dfaList=flucFunc(timeSeries,lag,qList,order)

    # Here I am extracting the Hurst Parameter for each fluctuation function:.
    # I use the transpose of dfaList for the purpose of easily parsing each row of the matrix
    h=hList(lag,dfaList.T, qList)

    # Now find and plot the Renyi exponent spectrum:
    tauFunc=renyiExp(qList,h)

    # Finally, from the Renyi exponents, find and plot the multifractal spectrum:
    mfList=mfSpec(tauFunc, qList)

    # Here we are calculating the PSD function:
    # Choose 22 Hz as the sampling frequency, as in Gomes et al.
    fs = 22 # Hz
    fList, PSDlist = PSDfunc(timeSeries,fs) 

    return None

#############################################################################

#Our first helper function will allow us to find the fluctuation function
#of our data, as described in the Gomes article.

#Here, "order" represents the order of the polynomial fitting,
# "qList" represents the list of powers of the fluctuation function, (q=2 for a monofractal case)
# "lag" represents the list of lags as a series of integers that represent side lengths
# bmag is the total magnitude of the magnetic field data

def flucFunc(bmag, lag, qList, order):
    
    # Here we use the MFDFA function to give us the fluctuation function ,
    # which we then take the log of:
    lag, dfaList = MFDFA(bmag, lag = lag, q = qList, order = order)
    
    #Obtain the log-log plot which has a slope of the log of the Hurst parameter:
    # This will plot the set of points for each lag s.
   
    #plt.loglog(lag,dfaList,'o')
    
        #label=f'fOU: MFDFA {qList[i]}=2'
    fig, ax = plt.subplots()
    for q in range(len(dfaList.T)):
        # Calculate the log of the DFA list and the lag:
        logdfa=np.log(dfaList.T[q])
        logs=np.log(lag)
        qArray=np.full(len(lag),qList[q])
        # Here we plot the fluctuation function with a color bar for each q value:
        sc = ax.scatter(logs,logdfa,c=qArray,s=8, cmap='inferno', vmin=-20, vmax=20)
    
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

# We also need a list of the Hurst parameters for varying q values:
def hList(lag,dfaList, qList):
    # Initialize an empty list to hold our Hurst parameters:
    hList=[]
    
    # For every q in out list of powers, we need to find the Hurst parameter:
    for dfa in dfaList:
        # We fit to a 1st degree polynomial for a linear fit
        H=np.polyfit(np.log(lag)[2:len(lag)],np.log(dfa[2:len(lag)]),1)[0]
        
        # print(np.polyfit(np.log(lag)[2:9],np.log(dfa[2:9]),1)[0])
        # print('Estimated H = '+'{:.3f}'.format(H))

        hList.append(H)
    
    plt.plot(qList,hList)
    plt.ylabel('h(q)')
    plt.xlabel('q')
    plt.title('Hurst Exponent Spectrum')
    plt.show()

    return hList


#############################################################################

# To find the Renyi exponents, we first need a list of Hurst parameters
# associated with q values.
def renyiExp(qList,hList):
    #Initialize an empty list to hold the Renyi exponents:
    tauList=[]
    
    # Now, for each q and H, we can calculate the Renyi exponent defined by
    # the above equation 
    for i in range(len(qList)):
        # Here we apply the relation of the Hurst parameter to the Renyi exponent:
        tauList.append(qList[i]*hList[i]-1)
    
    # Plot the resulting spectrum:
    plt.plot(qList,tauList,'o')
    plt.xlabel(r'$q$')
    plt.ylabel(r'$\tau(q)$')
    plt.title('Renyi Exponent Spectrum')
    plt.show()
    
    return tauList        

#############################################################################

def mfSpec(tauFunc,qList):
    
    # First we create a list of derivatives of the Renyi spectrum:
    alphaList=np.gradient(tauFunc,qList)

    #Initialize an empty list that will hold the f(\alpha) values:
    mfList=[]
    for i in range(len(tauFunc)):
        #Perform the operation written in markdown above this code block:
        mfList.append(qList[i]*alphaList[i]-tauFunc[i])

    # Finally, plot the f(\alpha) function
    plt.plot(alphaList,mfList)
    plt.title('Singularity Spectrum')
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$f(\alpha)$')
    plt.show()

    return mfList

#############################################################################

# This function takes in the bmag data and the sampling frequency.
def PSDfunc(bmag, fs):
    f, PSD = signal.periodogram(bmag,fs)
    plt.loglog(f,PSD)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power Spectral Density [nt**2/Hz]')  

    #print(len(f))
    #print((len(PSD)))
    
    # Now we take the slope of the log-log plot for varying frequency regimes.
    
   
    # UNCOMMENT THIS CODE TO FIT A PORTION OF THE DATA:
    """
    # For the fitting we also need to crop the PSD data over a frequency band:        
    freqBand1 = f[int(len(f)/500):int(len(f)/50)]
    PSDband1 = PSD[int(len(f)/500):int(len(f)/50)]

    # Find the linear fit of the selected portion of data, extract the slope and intercept:
    linearFit=np.polyfit(np.log(freqBand1),np.log(PSDband1),1)  
    slope=linearFit[0] 
    intercept=linearFit[1]

    # Plot the linear fit on a log-log scale
    plt.loglog(freqBand1, np.exp(intercept)*freqBand1**slope)
    plt.show()
    """


    return f,PSD

#############################################################################

# This function returns the shuffled data for analysis:
def shuffle(timeSeries):
    to_shuffle = np.copy(timeSeries)
    rng = np.random.default_rng()
    rng.shuffle(to_shuffle)
    #print(to_shuffle)

    return to_shuffle

#############################################################################

# This is simply a function to calculate the magnitudes of vectors in a list.
# The purpose is to handle PSP magnetic field data, which is only downloadable# in RTN (radial-tangential-normal) coordinates. We want the magnitude of the # magnetic field which is calculated as Bmag = sqrt(B_R**2+B_T**2+B_N**2)

def magnitude(list):

    bmagList = []

    for vector in list:
        # Calculate the magnitude here:
        bmag = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)  
        bmagList.append(bmag)
  
    return np.array(bmagList)

