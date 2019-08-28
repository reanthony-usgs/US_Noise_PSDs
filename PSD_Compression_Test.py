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

test_run = True

path = '/data/TA_ALASKA_2018_BHZ_PSDS/'

net, chan = "TA", "LHZ"

# Size of PSD Window in Seconds
window = 3600
Comp_Tres = 0.9

# This needs changed for different sample rates 
# Suggested BHZ - 2**15, LHZ - 2**10
nfft = 2**10
windlap = 0.75

##################################### Psd calculation
def get_dataless(net, chan, stime, etime, client):
    if os.path.exists(path + net + '_metadata.pickle'):
        with open(path + net + '_metadata.pickle', 'rb') as fhand:
            inv = pickle.load(fhand)
    else:
        client = Client()
        inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                                  channel=chan, network=net, level="response")
        with open(path + net + '_metadata.pickle', 'wb') as fhand:
            pickle.dump(inv, fhand)
    return inv



def write_results(net, sta, loc, chan, ctime, power, freq):
    if not os.path.exists(path + sta + '_PSD'):
        os.mkdir(path + sta + '_PSD')

    if not os.path.exists(path + sta + '_PSD/' + net + '_' + sta + '_' + loc + '_' + chan + '_freqs.txt'):
        f = open(path + sta + '_PSD/' + net + '_' + sta + '_' + chan + '_freqs.txt', 'w')
        for fr in freq:
            f.write(str(fr) + '\n')
        f.close()
    fname = 'PSD_' + net + '_' + sta + '_' + loc + '_' + chan + '_'
    fname += str(ctime.year) + '_' + str(ctime.julday).zfill(3)
    fname += str(ctime.hour).zfill(2)
    f = open(path + sta + '_PSD/' + fname+ '.pckl', 'wb')
    pickle.dump(power, f)
    f.close()
    return

def calc_psd(net, sta, chan, ctime, inv_sta, debug = False):
    estime = ctime + (24.*60.*60.)
    result_str =  str(ctime.julday) + ', ' + str(ctime.year) + ', ' + chan + '\n'
    try:
        st = client.get_waveforms(net, sta, "*", chan, ctime, estime,
                                   attach_response=False)
        st.detrend('constant')
        st.merge(fill_value=0.)
        # This should make sure we have 24 complete hours of data
        st.trim(starttime=ctime, endtime=estime, pad=True, fill_value=0.)
        if debug:
            print(st)
    except:
        return 'IRIS Problem, ' + result_str

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
                resp_inv = inv_sta.get_response(stT[0].id, stT[0].stats.starttime)
                resp, freqR = resp_inv.get_evalresp_response(t_samp=1./fs,
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


    return 'Complete, ' + result_str


# Grab station start time and end time
stime = UTCDateTime('1995-001T00:00:00')
etime = UTCDateTime('2019-001T00:00:00')

if test_run:
    client = Client()
    inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                              channel=chan, network=net, level="response")
    print(inv)
else:
    inv = get_dataless(net, chan, stime, etime, client)



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
        etime = stime + 365*24*60*60
    ctime = stime
    f = open(path + 'log_file_' + tsta, 'w')
    while ctime <= etime:
        try:
            info = calc_psd(tnet, tsta, chan, ctime, inv_sta)
            f.write(info)
        except:
            f.write('Outerloop issue, ' + str(ctime.julday) + ', ' + str(ctime.year) + ', ' + chan + '\n')
        ctime += 24.*60.*60.
    f.close()
    return

#for net_sta in nets_stas:
#    print('on sta: ' + net_sta.split('_')[1])
#    run_station(net_sta)
from multiprocessing import Pool
pool = Pool(30)
pool.map(run_station, nets_stas)
