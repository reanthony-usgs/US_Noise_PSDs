#!/usr/bin/env python
import glob
import pickle 
import numpy as np
from obspy.core import UTCDateTime

#### PSD Statistics ##############


def get_psd_stats(DB_path, sta, starttime, endtime, prctile, debug = True):
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
        pickle_file = open(path, 'rb')
        a = pickle.load(pickle_file)
        print(a)
        pickle_file.close()
        
        
        
        #directs = path.split('/')
        #f = directs[-1]
        #info = f.split('_')
        #if debug:
            #print(info)
        #year, time = info[5], info[6] 
        #day, hour = time[0:3],time[3:5]
        
        #PSDtime = UTCDateTime(year = int(year), julday = int(day), hour = int(hour))

        #if PSDtime >= starttime and PSDtime <= endtime:
            ##try:
            #if True:
                #pickle_file = open(path, 'rb')
                #a = pickle.load(pickle_file)
                
                #PSD_Array.append(a)
                #pickle_file.close()
            ##except:
            ##    print('Problem with:' + path)
    #print(PSD_Array)
    #PSD_Array = np.asarray(PSD_Array)    
    #prctile_psd = np.percentile(PSD_Array, prctile,0)
    
    #F_path = glob.glob(DB_path +  sta + '_PSD/*_freqs.txt')
    #freqs = np.loadtxt(str(F_path[0]))
    
    #return 1./freqs, prctile_psd
    return
