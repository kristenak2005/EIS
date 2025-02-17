# %% [markdown]
# **Python code for downloading data from VSO**

# %%
import os
from sunpy.net import Fido, attrs as a
from sunpy.time import parse_time
from astropy import units as u
import datetime as dt
from scipy.io import readsav

# times are changed for each file, one minute of each file. Example:
file_start_time = '2015/10/18 10:26:00'
file_end_time = '2015/10/18 10:27:00'
file_date = dt.datetime.strftime(dt.datetime.strptime(file_start_time,'%Y/%m/%d %H:%M:%S'), '%Y/%m/%d')
cadence = a.Sample(10*u.second)

# Where to put the data
if os.uname().sysname == 'Darwin':
    data_directory = '/mnt/scratch/data/spruksk2/SDO/'
else:
    data_directory = '/mnt/scratch/data/spruksk2/SDO/'

# Define date
date = parse_time(file_start_time)
attrs_time = a.Time(file_start_time, file_end_time)


# %% [markdown]
# For downloading AIA data

# %%
## AIA: passband changed to 193 A
instr = a.Instrument('AIA')
provider = a.Provider('JSOC')
passband = ['193']
## Define the output location & create the output directory if it doesn't exist
for pband in passband:
    print('Now processing '+pband+'A passband')
    data_dir = data_directory+pband.rjust(4, "0")+'/'+file_date
    os.makedirs(data_dir, exist_ok='True')
    wvlnth = a.Wavelength(int(pband)*u.Angstrom, int(pband)*u.Angstrom)
    result = Fido.search(attrs_time, instr, wvlnth, cadence)
    files = Fido.fetch(result, path = data_dir, overwrite=False, progress=True)


# %% [markdown]
# For downloading HMI data

# %%
## HMI:
instr = a.Instrument('HMI')
obs = a.Physobs.los_magnetic_field
#obs = a.Physobs.intensity
## Define the output location & create the output directory if it doesn't exist
data_dir = data_directory+'HMI_mag/'+file_date
#data_dir = data_directory+'HMI_cont/'+file_date
os.makedirs(data_dir, exist_ok='True')
result = Fido.search(attrs_time, instr, obs, cadence)
files = Fido.fetch(result, path=data_dir, overwrite=False, progress=True)



