#!/usr/bin/env python
import time
import pickle
import os
import glob
import numpy as np
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from scipy.signal import welch
from matplotlib.mlab import csd

# "Global" variables
client = Client("IRIS")

test_run = True

path = '/data/test_case/'

net, chan1, chan2 = "US", "BHZ", "BH1"

# Size of PSD Window in Seconds
window = 3600
Overlap = 0.5

# QC Parametres 
Comp_Tres = 0.9
SM_SF = 0.125
SM_EF = 0.25


# This needs changed for different sample rates !! 
# Suggested BHZ - 2**15, LHZ - 2**10
nfft = 2**15
windlap = 0.75

##################################### Psd calculation
def get_dataless(net, stime, etime, client):
    if os.path.exists(path + net + '_metadata.pickle'):
        with open(path + net + '_metadata.pickle', 'rb') as fhand:
            inv = pickle.load(fhand)
    else:
        client = Client()
        inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                                  channel="*H*", network=net, level="response")                    
        with open(path + net + '_metadata.pickle', 'wb') as fhand:
            pickle.dump(inv, fhand)
    return inv


def write_results(net, sta, loc, chan1, chan2, ctime, power, freq, phase=False):
    if phase:
        tag = 'PH'
    else:
        tag = 'PSD'
    if not os.path.exists(path + sta + '_' + tag):
        os.mkdir(path + sta + '_' + tag)

    if not os.path.exists(path + sta + '_'+ tag +'/' + net + '_' + sta + '_' + 
                          loc + '_' + chan1 + '_' + chan2+ '_freqs.txt'):
        f = open(path + sta + '_' + tag + '/' + net + '_' + 
                 sta + '_' + chan1 + '_' + chan2 +  '_freqs.txt', 'w')
        for fr in freq:
            f.write(str(fr) + '\n')
        f.close()
    fname = tag + '_' + net + '_' + sta + '_' + loc + '_' + chan1 + '_' + chan2 + '_'
    fname += str(ctime.year) + '_' + str(ctime.julday).zfill(3)
    fname += str(ctime.hour).zfill(2)
    f = open(path + sta + '_' + tag + '/' + fname+ '.pckl', 'wb')
    pickle.dump(power, f)
    f.close()
    return
    
    
def psd_done(net, sta, loc, chan1, chan2, ctime):
    fname = 'PSD_' + net + '_' + sta + '_' + loc + '_' + chan1 + '_' + chan2 + '_'
    fname += str(ctime.year) + '_' + str(ctime.julday).zfill(3)
    fname += str(ctime.hour).zfill(2)
    print(fname)
    return os.path.exists(path + sta + '_PSD/' + fname + '.pckl')
    
def data_done(net, sta, chan1, chan2, ctime):
    fname = 'PSD_' + net + '_' + sta + '_*_' + chan1 + '_' + chan2 + '_'
    fname += str(ctime.year) + '_' + str(ctime.julday).zfill(3)
    files = glob.glob(path + sta + '_PSD/' + fname + '*.pckl')
    return len(files)
        

def calc_psd(net, sta, chan1, chan2, ctime, inv_sta, debug = False):
    estime = ctime + (24.*60.*60.)
    result_str =  str(ctime.julday) + ', ' + str(ctime.year) + ', ' + chan1 + ', ' + chan2 + '\n'
    try:

        num_files = data_done(net, sta, chan1, chan2, ctime)
        if num_files == 24:
            return 'Not grabbing'
        st = client.get_waveforms(net, sta, "*", chan1, ctime, estime,
                                   attach_response=False)
        
        st += client.get_waveforms(net, sta, "*", chan2, ctime, estime,
                                   attach_response=False)
        
        st.detrend('constant')
        st.merge(fill_value=0.)
        # This should make sure we have 24 complete hours of data
        st.trim(starttime=ctime, endtime=estime, pad=True, fill_value=0.)
        if debug:
            print(st)
    except:
        return 'IRIS Problem, ' + result_str
    
    for stT in st.slide(window, window*overlap):
        if psd_done(net, sta, stT[0].stats.location, chan1, chan2, stT[0].stats.starttime):
            if debug:
                print('Skipping:' + net + ' ' + sta + ' ' + stT[0].stats.location +  ' ' + chan1, + ' ' + chan2)
                print(ctime)
            continue
        if 'fs' not in vars():
                fs = 1./stT[0].stats.delta
        if chan1 == chan2:
            data_count = np.count_nonzero(stT[0].data)
        else:    
            data_count1 = np.count_nonzero(stT[0].data)
            data_count2 = np.count_nonzero(stT[1].data)
            data_count = min([data_count1, data_count2])
        total_count = window*fs
        if data_count/total_count >= Comp_Tres:
            try:
                # Default uses a Hann Window and demeans the data
                # Need to switch to cross-power
                if chan1 == chan2:
                    idx1, idx2 = 0, 0
                else:
                    idx1, idx2 = 0, 1
                cxy, freq = csd(stT[idx1].data, stT[idx2].data, NFFT=nfft, 
                                noverlap = nfft*windlap, Fs=fs,
                                scale_by_freq=True, detrend='linear')
                
                # Be careful with cxy it is complex and not real
                freq, cxy = freq[1:], cxy[1:]
                # if this is a time sink we could change this
                resp_inv1 = inv_sta.get_response(stT[0].id, stT[0].stats.starttime)
                loc = stT[0].stats.location
                if chan1 == chan2:
                    resp, freqR = resp_inv1.get_evalresp_response(t_samp=1./fs,
                                                             nfft=nfft, output='ACC')
                    resp = resp[1:]
                    power = np.abs(cxy)
                    power = 10.*np.log10(power/(np.abs(resp)**2))
                    
                else:
                    resp_inv2 = inv_sta.get_response(stT[1].id, stT[1].stats.starttime)
                    resp1, freqR = resp_inv1.get_evalresp_response(t_samp=1./fs,
                                                             nfft=nfft, output='ACC')
                    resp2, freqR = resp_inv2.get_evalresp_response(t_samp=1./fs,
                                                             nfft=nfft, output='ACC')
                    resp1, resp2 = resp1[1:], resp2[1:]
                    cxy /= ((np.abs(resp1)**2)*(np.abs(resp2)**2))
                    cxy *=  resp2*np.conjugate(resp1)
                    # Now we have the response removed
                    power = 10.*np.log10(np.abs(cxy))
                    phase = np.unwrap(np.angle(cxy))
                    write_results(net, sta, loc, chan1, chan2, stT[0].stats.starttime, power, freq, phase=True)
                write_results(net, sta, loc, chan1, chan2, stT[0].stats.starttime, power, freq)
                    
                
            except:
                print('Problem computing PSD' + stT[0].id())
                continue
        else:
            print('Not enough data in this window:', (data_count/total_count)*100,
            '% Complete. Skipping this Window.')
            continue


    return 'Complete, ' + result_str


# Grab station start time and end time
stime = UTCDateTime('1995-001T00:00:00')
etime = UTCDateTime('2019-001T00:00:00')

if test_run:
    client = Client()
    inv = client.get_stations(starttime=stime, endtime=etime, station="*K",
                              channel="*H*", network=net, level="response")
    print(inv)
else:
    inv = get_dataless(net, stime, etime, client)



nets_stas = []
for net in inv:
    for sta in net:
        nets_stas.append(str(net.code + '_' + sta.code))


print(nets_stas)


def run_station(net_sta):
    tnet, tsta = net_sta.split('_')
    # Using inv as a global variable
    inv_sta = inv.select(station=tsta)
    stime = inv_sta[0][0].start_date
    etime = inv_sta[0][0].end_date
    if etime >= UTCDateTime('2019-001T00:00:00'):
        etime = UTCDateTime('2019-001T00:00:00')
    if test_run:
        stime = UTCDateTime('2018-001T00:00:00')
        etime = stime + 10*24*60*60
    ctime = stime
    if not os.path.exists(path):
        os.mkdir(path)
    f = open(path + 'log_file_' + tsta, 'w')
    while ctime <= etime:
        try:
            info = calc_psd(tnet, tsta, chan1, chan2, ctime, inv_sta)
            if info == 'Not grabbing':
                pass
            else:
                f.write(info)
        except:
            f.write('Outerloop issue, ' + str(ctime.julday) + ', ' + str(ctime.year) + ', ' + chan1 + ', '  + chan2 + '\n')
        ctime += 24.*60.*60.
    f.close()
    return

#for net_sta in nets_stas:
#    print('on sta: ' + net_sta.split('_')[1])
#    run_station(net_sta)
from multiprocessing import Pool
pool = Pool(10)
pool.map(run_station, nets_stas)
