#!/usr/bin/env python

import time
import pickle
import os

import numpy as np

from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from scipy.signal import welch

# "Global" variables
client = Client("IRIS")
debug = True
time_run = True

# Size of PSD Window in Seconds
window = 3600
Comp_Tres = 0.9

nfft = 2**15
windlap = 0.75

net, sta, chan = "US", "ECSD", "BHZ"

##################################### Psd calculation


def write_results(net, sta, loc, chan, ctime, power, freq):
    if not os.path.exists(net + '_' + sta + '_' + loc + '_' + chan + '_freqs.txt'):
        f = open(net + '_' + sta + '_' + chan + '_freqs.txt', 'w')
        for fr in freq:
            f.write(str(fr) + '\n')
        f.close()
    fname = 'PSD_' + net + '_' + sta + '_' + loc + '_' + chan + '_'
    fname += str(ctime.year) + '_' + str(ctime.julday).zfill(3)
    fname += str(ctime.hour).zfill(2)
    if debug:
        fname += str(ctime.minute).zfill(2)
    f = open(fname+ '.pckl', 'wb')
    pickle.dump(power, f)
    f.close()
    return


def calc_psd(net, sta, chan, ctime, time_run):
    if time:
        data_time = time.time()
    estime = ctime + (24.*60.*60.) 
    
    try:
        st = client.get_waveforms(net, sta, "*", chan, ctime, estime,
                                   attach_response=True)
        st.detrend('constant')                           
        st.merge(fill_value=0.)    
                    
        # This should make sure we have 24 complete hours of data
        st.trim(starttime=ctime, endtime=estime, pad=True, fill_value=0.)
        if debug:
            print(st) 
        
    except:
        print('IRIS issue for: ' + net + ' ' + sta + str(ctime))
        return
    if time:
        print('It took', time.time() - data_time , ' s to load in all the data')
        psd_time = time.time()


    for stT in st.slide(window, window):
        
        if 'fs' not in vars():
                fs = 1./stT[0].stats.delta
        
        data_count = np.count_nonzero(stT[0].data)
        total_count = window*fs
        
        if data_count/total_count >= Comp_Tres:
            try:
                # Default uses a Hann Window and demeans the data
                freq, power = welch(stT[0].data, fs=fs, nperseg=nfft,
                                    noverlap=nfft*windlap, detrend='linear')
                freq, power = freq[1:], power[1:]
                # if this is a time sink we could change this
    
                resp, freqR = stT[0].stats.response.get_evalresp_response(t_samp=1./fs,
                                                                      nfft=nfft, output='ACC')
                resp = resp[1:]
                power = np.abs(power)
                power = 10.*np.log10(power/(np.abs(resp)**2))
                loc = stT[0].stats.location
                write_results(net, sta, loc, chan, stT[0].stats.starttime, power, freq)
                
                
            except:
                print('Problem computing PSD' + stT[0].id())
                continue
        else:
            print('Not enough data in this window:', (data_count/total_count)*100,
            '% Complete. Skipping this Window.') 
            continue
            
    if time:
        print('It took', time.time() - psd_time , 's to calculate the PSDs')
        
    return
        

# Grab station start time and end time
stime = UTCDateTime('1995-001T00:00:00')
etime = UTCDateTime('2019-001T00:00:00')
inv = client.get_stations(starttime=stime, endtime=etime, station=sta,
                          channel=chan, network=net, level="response")

if debug:
    print(inv)
for netinv in inv:
    for stainv in netinv:
        stime = stainv.start_date
        etime = stainv.end_date

# Keep this in for testing and remove for actual run
stime = UTCDateTime('2017-001T00:00:00')
etime = UTCDateTime('2017-006T00:00:00')
# Loop over the times for each PSD calculation
ctime = stime

if time_run:
    total_time = time.time()

while ctime <= etime:
    calc_psd(net, sta, chan, ctime, time_run)
    ctime += 24.*60.*60.

if time_run:
    print('It took', time.time() - total_time , ' s to do the entire computation')
