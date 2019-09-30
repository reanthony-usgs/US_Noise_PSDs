#!/usr/bin/env python

# Ran a test run on compute_psd.py for ANMO requesting the data from IRIS. 

# For Verification Purposes, rerun the same data locally and calulate a PSD

# This script reruns the PSD and plots the output as well as the PSD from the database


from obspy.core import Stream, read, UTCDateTime, inventory 
from obspy.core.inventory.inventory import read_inventory 
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np
from scipy.signal import welch
from obspy.signal.invsim import evalresp
from math import pi
import sys
from itertools import combinations
import math, copy
import random
import csv
import pickle 


debug = True

net = "IU"
sta = "ANMO"
loc = "00"
chan = "BHZ"

SDay = "005";
SHour = "18";
SMin = "30"

Path = "/data/test_case/ANMO_PSD/"


# Parameters for Code ##################################################

# input Sample Rate
fs = 40

# subwindow size and overlap
nfft = 2**15
windlap = 0.75 
# size of window and overlap for PSD
window = 3600
overlap = 0.5

# Data completeness for each window
Comp_Thres = 0.9

########################################################################

seed_ID = net + "." + sta + "." + loc + "." + chan
TestHour = UTCDateTime("2019-" + SDay + "T" + SHour + ":" + SMin + ":00")


# Get the response

resppath = "/APPS/metadata/RESPS/RESP." + net + "." + sta + "." + loc + "." + chan
inv = read_inventory(resppath)
inv_Resp = inv.get_response(seed_ID,TestHour)
resp, freqR = inv_Resp.get_evalresp_response(t_samp=1./fs, nfft=nfft, output='ACC')
resp = resp[1:]

# read the first 11 days of ANMO data
st = Stream()

for day in range(1,12):
    try:
        st += read("/msd/" + net + "_" + sta + "/" + str(TestHour.year) + "/" + str(day).zfill(3) + "/" + loc + "_" + chan + "*.seed")
    except:
        print("No Data Here")
    
    
total_count = window*fs

st.trim(TestHour, TestHour+window)
if debug:
    st.plot()
st.detrend('constant')

st.merge(fill_value=0.)
data_count = np.count_nonzero(st[0].data)
        
if data_count/total_count >= Comp_Thres:
    
    tr = st[0]
    # Default uses a Hann Window and demeans the data
    freq, power = welch(tr.data, fs=fs, nperseg = nfft, 
                            noverlap = nfft*windlap, detrend = 'linear')
            
    freq = freq [1:]
    power = power[1:] 
            
    power = np.abs(power)
    power = 10.*np.log10(power/(np.abs(resp)**2))
            
    # get rid of the extra digits
    power = np.round(power,10)
else:
    print('Not Enough Data')    
            

# Load and unpickle the compute_psd.py PSD for the same Hour 

PSDpath = (Path + 'PSD_' + net + '_' + sta + '_' + loc + "_BHZ_BHZ_2019_" + SDay + SHour + SMin + '.pckl' )
Fpath = (Path + 'IU_ANMO_BHZ_BHZ_freqs.txt')

freq2 = np.loadtxt(Fpath)


pickle_file = open(PSDpath, 'rb')
power2 = pickle.load(pickle_file)
 
if debug:
    print(np.shape(power2))
    print(PSDpath)

# Plot up the PSDs

NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()


fig = plt.figure(2)

plt.semilogx(1./freq2, power2, linewidth=5., color='r', label='TA Algorithm')
plt.semilogx(1./freq, power, linewidth=5., color='b', label='Local Algorithm')
plt.semilogx(NLNMper, NLNMpower, linewidth=2., color='k')
plt.semilogx(NHNMper, NHNMpower, linewidth=2., color='k',label='NLNM/NHNM')
plt.xlabel('Period (s)')
plt.ylabel('Power (dB)')
plt.xlim((0.05,500.))
plt.ylim((-200., -80.))
plt.legend()
plt.show()                     
