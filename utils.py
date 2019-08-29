#!/usr/bin/env python
import glob
import pickle
import os
import numpy as np
from obspy.core import UTCDateTime

#### PSD Statistics ##############

path = '/data/LHZ_TA/'

def get_psd_stats(DB_path, sta, starttime, endtime, prctile, debug = False):
    """
    Inputs:
    DB_path: Path to the directory 
    sta: station name
    starttime and endtime: start and end times
    prctile: the percentile you want 

    Outputs
    Periods, Power in dB rel. 1 (m/s^2)^2/Hz.
    """
    
    starttime, endtime = UTCDateTime(starttime), UTCDateTime(endtime)
    All_PSD_Paths = glob.glob(DB_path + sta + '_PSD/*.pckl') 
    if debug:
        print(All_PSD_Paths)
    PSD_Array = []
    for path in All_PSD_Paths:
        directs = path.split('/')
        f = directs[-1]
        info = f.split('_')
        if debug:
            print(info)
        year, time = info[5], info[6] 
        day, hour = time[0:3],time[3:5]
        PSDtime = UTCDateTime(year = int(year), julday = int(day), hour = int(hour))
        if PSDtime >= starttime and PSDtime <= endtime:
            pickle_file = open(path, 'rb')
            a = pickle.load(pickle_file)
            PSD_Array.append(a)
            pickle_file.close()
    if debug:
        print(PSD_Array)
    PSD_Array = np.asarray(PSD_Array)    
    prctile_psd = np.percentile(PSD_Array, prctile,0)    
    F_path = glob.glob(DB_path +  sta + '_PSD/*_freqs.txt')
    freqs = np.loadtxt(str(F_path[0]))
    pers = np.array(1./freqs)
    prctile_psd = np.array(prctile_psd)
    return 1./freqs, prctile_psd
    
def get_dataless(net, chan, stime, etime, client):
    if os.path.exists(path + net + '_metadata.pickle'):
        with open(path + net + '_metadata.pickle', 'rb') as fhand:
            inv = pickle.load(fhand)
    else:
        
        inv = client.get_stations(starttime=stime, endtime=etime, station="*K",
                                  channel=chan, network=net, level="response")
        with open(path + net + '_metadata.pickle', 'wb') as fhand:
            pickle.dump(inv, fhand)
    return inv
