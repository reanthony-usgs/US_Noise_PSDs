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

path = '/home/aringler/data_stuff/US_array/'
sta = 'J20K'
stime = UTCDateTime('2019-161T00:00:00')
etime = UTCDateTime('2019-172T00:00:00')
per, psd = utils.get_psd_stats(path, sta, stime, etime, 10.)







#def setupmap(central_lon, central_lat,handle):
    ##handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    #handle.set_extent(extent)

    #handle.add_feature(cfeature.LAND)
    #handle.add_feature(cfeature.OCEAN)
    #handle.add_feature(cfeature.COASTLINE)
    #handle.add_feature(cfeature.BORDERS, linestyle=':')
    #handle.add_feature(cfeature.LAKES, alpha=0.5)
    #handle.add_feature(cfeature.RIVERS)
    #handle.add_feature(cfeature.STATES, edgecolor='gray')
    #return handle

