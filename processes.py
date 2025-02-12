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

fitted_lines = {
    "fe_12_195" : ["fe_12_195_119.2c.template.h5",0,"Fe XII 195$\\AA$"],
    "ar_14_194" : ["ar_14_194_396.6c.template.h5",5,"Ar XIV 194$\\AA$"],
    "ca_14_193" : ["ca_14_193_874.6c.template.h5",1,"Ca XIV 193$\\AA$"],
    "si_10_258" : ["si_10_258_375.1c.template.h5",0,"Si X 258$\\AA$"],
    "s_10_264" : ["s__10_264_233.1c.template.h5",0,"S X 264$\\AA$"],
    "fe_13_202" : ["fe_13_202_044.1c.template.h5",0,"Fe XIII 202$\\AA$"],
    "fe_13_203" : ["fe_13_203_826.2c.template.h5",1,"Fe XIII 203$\\AA$"]
}


## Function to get the list of AIA filenames
def get_aia_filelist(data_dir, passband, file_date):
    files = glob.glob(data_dir+str(passband).rjust(4, "0")+'/'+file_date+'/*.fits', recursive=True)
    files.sort()
    files_dt = []
    for file_i in files:
        hdr = fits.getheader(file_i, 1)
        try:
            files_dt.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%fZ'))
        except:
            files_dt.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%f'))
    return files, files_dt




