#!/usr/bin/env python

# Purpose of this code is to:

# 1) time how long it takes a single compute node
# to request and process 1-year of TA data on a single channel 

# 2) Estimate the size of a full TA PSD dataset without compression

# 3) Estimate the size of a full TA PSD dataset in a python Pickle 


# We will use Welch's Method of section averaging on 1-hour windows of data using
# 13-minute subsections with 75% overlap. 

# This will be similar to peterson's algorithm in that 15 subsections will be used

# Hann windows and demeaning and detrending will be performed on each subsegment. 
# The "linear" detrend option in SciPy Welch accomplishes this. This subsegment
# Processing is identical to Petersons (1993) method


# This should be a nice tradeoff between resolution out to periods of interest 
# for these stations (a few hundred s period) while not being overly dominated
# by teleseisms in the primary microseism band

# This parameters are similar to Petersons, so comparison of the PSDs
# with the NHNM and NLNM seems appropriate 


from obspy.core import read, UTCDateTime
from scipy.signal import welch
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
from obspy.signal.invsim import evalresp
import csv
import pickle 

from obspy.clients.fdsn import Client
client = Client("IRIS")


debug = True

# Do you want timing information spit out?
Time = True


# Set all our parameters 

# Size of PSD Window in Seconds 
window = 3600

nfft = 2**15
windlap = 0.75

# start our timer
start_time = time.time()

# Load in data 
st = read()
st.clear()

for day in range(5,106):
    
    print(day)

    # Set start and end times 
    
    if day == 365:
        stime = UTCDateTime('2017-365T00:00:00.000')
        etime = UTCDateTime ('2018-001T00:00:00.000')
        
    else:
        stime = UTCDateTime('2019-' + str(day).zfill(3) + 'T00:00:00.000')
        etime = UTCDateTime('2019-' + str(day+1).zfill(3) + 'T00:00:00.000')
    
    try: 
        st += client.get_waveforms("IU", "ANMO", "00", "BHZ", stime, etime, attach_response=True)
        
    except:
        print('Unknown IRIS Failure')
        
data_time = time.time()     
                
if Time:
    print('It took', data_time - start_time, ' s to load in all the data')
    print(st)


# Calculate PSDs 

for window in st.slide(window,window):
    
    # Get all timing information 
    
    if debug:
        print(window[0].stats.starttime)
        
    try:
        
        tr = window[0]
        fs = 1./tr.stats.delta
        
        # Get all timing information 
        
        year = tr.stats.starttime.year
        julday = tr.stats.starttime.julday
        hour = tr.stats.starttime.hour
        minute = tr.stats.starttime.minute
        second = tr.stats.starttime.second
        
        # Default uses a Hann Window and demeans the data
        freq, power = welch(tr.data, fs=fs, nperseg = nfft, 
                        noverlap = nfft*windlap, detrend = 'linear')
        
        freq = freq [1:]
        power = power[1:] 
        
        power = np.abs(power)
        
        if debug:
            print('The size of the PSD is',np.shape(power))
        
        # Try to get the response which we attached 
        
        resp, freqR = tr.stats.response.get_evalresp_response(t_samp=tr.stats.delta,
                nfft=nfft, output='ACC')
        
        resp = resp[1:]
        
        if debug:
            print('Here is our Response:',resp)
            
        # Convert to dB and remove the response 
        
        # This has the added benefit of failing and throwing an except
        # of the data segment is not long enough
    
        Power = 10.*np.log10(power/(np.abs(resp)**2))
        
        # Assemble timing information Vector
        date_vec = [year, julday, hour, minute, second]
    
        # Assemble the Matrix of PSDs and dates
    
    
        if "Power_M" not in vars():
            Power_M = Power
            date_vec_M = date_vec
            
        else:
            Power_M = np.vstack((Power_M,Power))
            date_vec_M = np.vstack((date_vec_M,date_vec))
    except:
        print('Not Enough Data')


psd_time = time.time()  

if Time:
    print('It took', psd_time-data_time, 's to calculate the PSDs')
        
        
# Save power and date information as both a flat file and as a Pickle 
 
# Flat File 
freq = np.matrix.transpose(freq)

f=open('BHZ_freqs_ANMO.txt','w')
for fw in freq:
    f.write(str(fw) + '\n')
f.close()

f2=open('PSD_Times_ANMO.txt','w')
writer = csv.writer(f2, delimiter='\t')
writer.writerows(date_vec_M)
f2.close()

f3=open('PSDs_ANMO.txt','w')
writer = csv.writer(f3, delimiter='\t')
writer.writerows(Power_M)
f3.close()

flat_time = time.time()

if Time:
    print('It took', flat_time-psd_time, 's to write PSDs as flat files')
   
# pickle
filename = 'BHZ_freqs_ANMO.pckl'
outfile = open(filename, 'wb')
pickle.dump(freq, outfile)
outfile.close()

filename = 'PSD_Times_ANMO.pckl'
outfile = open(filename, 'wb')
pickle.dump(date_vec_M, outfile)
outfile.close()


filename = 'PSDs_ANMO.pckl'
outfile = open(filename, 'wb')
pickle.dump(Power_M, outfile)
outfile.close()

pickle_time = time.time()

if Time:
    print('It took', pickle_time-flat_time, 's to pickle PSDs')

                
    
    












