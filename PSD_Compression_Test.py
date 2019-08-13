#!/usr/bin/env python

import time
import pickle
import os

import numpy as np

from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from scipy.signal import welch

# Purpose of this code is to:

# 1) time how long it takes a single compute node
# to request and process 1-year of TA data on a single channel

# 2) Estimate the size of a full TA PSD dataset without compression

# 3) Estimate the size of a full TA PSD dataset in a python Pickle


# We will use Welch's Method of section averaging on
# 1-hour windows of data using 13-minute subsections with 75% overlap.

# This will be similar to peterson's algorithm in that
# 15 subsections will be used

# Hann windows and demeaning and detrending will be performed on
#each subsegment. The "linear" detrend option in SciPy Welch accomplishes
#this. This subsegment Processing is identical to Petersons (1993) method

# This should be a nice tradeoff between resolution out to periods of interest
# for these stations (a few hundred s period) while not being overly dominated
# by teleseisms in the primary microseism band

# This parameters are similar to Petersons, so comparison of the PSDs
# with the NHNM and NLNM seems appropriate

# "Global" variables
client = Client("IRIS")
debug = True

# Size of PSD Window in Seconds
window = 3600

nfft = 2**15
windlap = 0.75
Time = True

# start our timer
start_time = time.time()
net, sta, chan = "IU", "ANMO", "BHZ"

# error counters 

response_errors = 0
data_errors = 0


##################################### Psd calculation


def write_results(net, sta, chan, ctime, power, freq):
    # make a file of frequencies if they don't already exist 
    # If they do exist, skip
    if not os.path.exists(net + '_' + sta + '_' + chan + '_freqs.txt'):
        f = open(net + '_' + sta + '_' + chan + '_freqs.txt', 'w')
        for fr in freq:
            f.write(str(fr) + '\n')
        f.close()
    # write the psd file name     
    fname = 'PSD_' + net + '_' + sta + '_' + chan + '_'
    fname += str(ctime.year) + '_' + str(ctime.julday).zfill(3)
    fname += str(ctime.hour).zfill(2)
    f = open(fname + '.txt', 'w')
    for pow in power:
        f.write(str(pow) + '\n')
    f.close()
    f = open(fname+ '.pckl', 'wb')
    pickle.dump(power, f)
    f.close()
    return


def calc_psd(net, sta, chan, ctime, start_time, Time, data_errors, response_errors):
    start_time_data = time.time()
    #if True:
    try:
        st = client.get_waveforms(net, sta, "*", chan, ctime,
                                  ctime + (60.*60.), attach_response=True)
        st.detrend('constant')
        st.merge(fill_value=0.)
        data_time = time.time()
        if Time:
            print('It took', data_time - start_time_data, ' s to load in 1 Hour of data')
            print(st)
    except:
        print('IRIS issue for: ' + net + ' ' + sta + str(ctime))
        if debug:
            data_errors = data_errors + 1
        return 

    if 'fs' not in vars():
        fs = 1./st[0].stats.delta
    # Default uses a Hann Window and demeans the data
    freq, power = welch(st[0].data, fs=fs, nperseg=nfft,
                        noverlap=nfft*windlap, detrend='linear')
    freq, power = freq[1:], power[1:]
    # if this is a time sink we could change this
    # This threw an error on ANMO 2019-01-23 T00:00:00
    
    try:
        resp, freqR = st[0].stats.response.get_evalresp_response(t_samp=1./fs,
                                                          nfft=nfft, output='ACC')
        resp = resp[1:]
        power = np.abs(power)
        power = 10.*np.log10(power/(np.abs(resp)**2))

        psd_time = time.time()
        
    except:
        print('Error removing response')
        
        if debug:
            response_errors = response_errors + 1
            
        return    
        
    if Time:
        print('It took', psd_time-data_time, 's to calculate a single PSD')

    write_results(net, sta, chan, ctime, power, freq)

    pickle_time = time.time()
    if Time:
        print('It took', psd_time-pickle_time, 's to write a single PSD to both pickle and text files')
    return

# Grab station start time and end time
stime = UTCDateTime('2019-005T00:00:00')
etime = UTCDateTime('2019-106T00:00:00')
inv = client.get_stations(starttime=stime, endtime=etime, station=sta,
                          channel="BH*", network=net, level="response")

if debug:
    print(inv)

    
# What if the start of the station is not exactly on an hour?     
# We don't want this unless we're going for broke. If overrides the starttime to
# The station starttime (e.g. - ANMO 6/30/2008 when we want 1/5/2019)
#for netinv in inv:
#    for stainv in netinv:
#        stime = stainv.start_date
#        etime = stainv.end_date
        
#if debug:
#    print('Here is our start time', stime)        
        

# Loop over the times for each PSD calculation
ctime = stime
while ctime <= etime:
    if debug:
        print(ctime)
    calc_psd(net, sta, chan, ctime, start_time, Time, data_errors, response_errors)
    ctime += 60.*60.

# get a final time to do this

End_time = time.time()

if time:
    print('Total Time to run is', End_time-start_time, 's')
    
    
if debug:
    print('Total Number of errors reading in Data:', data_errors)
    print('Total Number of errors Removing Response:', response_errors)    
