#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from obspy.core import UTCDateTime
import utils
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
from obspy.clients.fdsn import Client

path = '/data/LHZ_TA/'
sta = 'A19K'
net ='TA'
chan ='LHZ'
client = Client('IRIS')
stime = UTCDateTime('2018-150T00:00:00')
etime = UTCDateTime('2018-151T00:00:00')
per, psd = utils.get_psd_stats(path, sta, stime, etime, 10.)



def setupmap(central_lon, central_lat,handle):
    #handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    handle.set_extent(extent)

    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle

inv = utils.get_dataless(net, chan, stime, etime, client)
lats, lons = [], []
for net in inv:
    for sta in net:
        for chans in sta:
            print(chans)
            #coors = inv.get_coordinates(chans)
            lats.append(chans.latitude)
            lons.append(chans.longitude)
lats = np.array(lats)
lons = np.array(lons)


fig= plt.figure(figsize=(12,12))


boxcoords=[min(lats) -1., min(lons)-1., max(lats) +1. , max(lons) + 1.]
extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])            
ax = plt.subplot(1,1,1, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, transform=ccrs.PlateCarree())
plt.show()
