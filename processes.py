import datetime as dt
import os
import matplotlib.pyplot as plt
from sunpy.map import Map
import numpy as np
import matplotlib.colors as colors
import eispac.net
from eispac.net.attrs import FileType
import glob
from astropy import units as u
from sunpy.time import parse_time
from astropy.io import fits as fits
from sunpy.net import Fido, attrs as a
from aiapy.calibrate import register, update_pointing, correct_degradation
import matplotlib.gridspec as gridspec
import sunkit_image.enhance as enhance
from astropy.coordinates import SkyCoord
from general_routines import closest
import asdf

#Specific timestamps for analysis
eis_evts = ['20151018_102719', '20151018_173743', '20151018_113839', '20151018_191443', '20151018_124939']


time_range = ['2015/10/18T10:27:19', '2015/10/18T12:49:39']
date = '20151018'

timerange = [dt.datetime.strptime(tr,'%Y/%m/%dT%H:%M:%S') for tr in time_range]
file_date = dt.datetime.strftime(dt.datetime.strptime(time_range[0],'%Y/%m/%dT%H:%M:%S'), '%Y/%m/%d')

