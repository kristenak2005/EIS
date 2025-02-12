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

#data_location contains h5 files
#+file_date is removed not needed
data_location = '/mnt/scratch/data/spruksk2/data'
#data_location = '/Users/dml/Data/EIS/'+file_date
#os... commented out as directory already exists
#os.makedirs(data_location, exist_ok=True)
#output_location = '/Users/dml/python_output/EIS_work/'+date
output_location = '/mnt/scratch/data/spruksk2/python_output/EIS_work/'+date
os.makedirs(output_location, exist_ok=True)

#results = Fido.search(a.Time(timerange[0], timerange[1]),
                    #  a.Instrument('EIS'),
                    #  a.Physobs.intensity,
                    #  a.Source('Hinode'),
                    #  a.Provider('NRL'),
                    #  a.Level('1'))  
#files = Fido.fetch(results, path = data_location, overwrite=False, progress=True)
files = sorted(glob.glob(os.path.join(data_location, '*.h5')))
#ar not included as not in data,
fitted_lines = {
    "fe_12_195" : [
        "eis_20151018_102719.fe_12_195_119.2c-0.fit.h5",
        "eis_20151018_124939.fe_12_195_119.2c-0.fit.h5",
        "eis_20151018_191443.fe_12_195_119.2c-0.fit.h5",
        "eis_20151018_113839.fe_12_195_119.2c-0.fit.h5",
        "eis_20151018_173743.fe_12_195_119.2c-0.fit.h5"
  ],
 "si_10_258" : [
    "eis_20151018_113839.si_10_258_375.1c-0.fit.h5",
    "eis_20151018_102719.si_10_258_375.1c-0.fit.h5",
    "eis_20151018_173743.si_10_258_375.1c-0.fit.h5",
    "eis_20151018_124939.si_10_258_375.1c-0.fit.h5", 
    "eis_20151018_191443.si_10_258_375.1c-0.fit.h5"
  ],
  "s_10_264" : [
    "eis_20151018_124939.s__10_264_233.1c-0.fit.h5",
    "eis_20151018_102719.s__10_264_233.1c-0.fit.h5",
    "eis_20151018_191443.s__10_264_233.1c-0.fit.h5",
    "eis_20151018_113839.s__10_264_233.1c-0.fit.h5",
    "eis_20151018_173743.s__10_264_233.1c-0.fit.h5"
  ],
  "fe_13_202" : [
    "eis_20151018_124939.fe_13_202_044.1c-0.fit.h5",
    "eis_20151018_102719.fe_13_202_044.1c-0.fit.h5",
    "eis_20151018_191443.fe_13_202_044.1c-0.fit.h5",
    "eis_20151018_113839.fe_13_202_044.1c-0.fit.h5",
    "eis_20151018_173743.fe_13_202_044.1c-0.fit.h5"
  ]
    
}




