# %% [markdown]
# **Notebook for downloading EIS data using EISPAC and producing a first look fitted image with an AIA & HMI context image**
# 
# First, we define the time period of interest and create the data and output folders.

# %%
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

# %%
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

# %%
## Function to get the list of HMI filenames
def get_hmi_filelist(data_dir, obs, file_date):
    files = glob.glob(data_dir+'HMI_mag/'+file_date+'/*.fits', recursive=True)
    files.sort()
    files_dt = []
    for file_i in files:
        hdr = fits.getheader(file_i, 1)
        try:
            files_dt.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%fZ'))
        except:
            files_dt.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%f'))
    return files, files_dt

# %% [markdown]
# Define function to plot the context images

# %%
def get_context_img(img_time):
# Define some constants.
    SDO_data_directory = '/mnt/scratch/data/spruksk2/SDO/'

    aia_passband = '193'
    cadence = a.Sample(12*u.second)
    aia_instr = a.Instrument('AIA')
    aia_provider = a.Provider('JSOC')
    aia_wvlnth = a.Wavelength(int(aia_passband)*u.Angstrom, int(aia_passband)*u.Angstrom)

    start_time = (parse_time(img_time)-10*u.second).strftime('%Y/%m/%d %H:%M:%S')
    end_time = (parse_time(img_time)+10*u.second).strftime('%Y/%m/%d %H:%M:%S')
    attrs_time = a.Time(start_time, end_time)
    file_date = (parse_time(img_time)-10*u.second).strftime('%Y/%m/%d')

# Download the AIA data    
    f_aia, file_time = get_aia_filelist(SDO_data_directory, aia_passband, file_date)
    if f_aia == []:
        print('Downloading data for '+str(aia_passband).rjust(4, "0")+' passband')
        data_dir = SDO_data_directory+aia_passband.rjust(4, "0")+'/'+file_date
        os.makedirs(data_dir, exist_ok='True')
        result = Fido.search(attrs_time, aia_instr, aia_wvlnth, cadence, aia_provider)
        files = Fido.fetch(result, path = data_dir, overwrite=False, progress=True)

    f_aia, file_time = get_aia_filelist(SDO_data_directory, aia_passband, file_date)
    nearest_img = closest(np.array(file_time), dt.datetime.strptime(start_time,"%Y/%m/%d %H:%M:%S"))
    if dt.datetime.strptime(start_time, "%Y/%m/%d %H:%M:%S") < file_time[nearest_img] < dt.datetime.strptime(end_time, "%Y/%m/%d %H:%M:%S"):
        print('Data already downloaded for '+str(aia_passband).rjust(4, "0")+' passband')
    else:
        print('Downloading data for '+str(aia_passband).rjust(4, "0")+' passband')
        data_dir = SDO_data_directory+aia_passband.rjust(4, "0")+'/'+file_date
        os.makedirs(data_dir, exist_ok='True')
        result = Fido.search(attrs_time, aia_instr, aia_wvlnth, cadence, aia_provider)
        files = Fido.fetch(result, path = data_dir, overwrite=False, progress=True)
        f_aia = glob.glob(data_dir+'/*.fits')
        f_aia.sort()

# Now find the nearest AIA image in time to the EIS map and load it as a map
    time = []
    for file_i in f_aia:
        hdr = fits.getheader(file_i, 1)
        try:
            time.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%fZ'))
        except:
            time.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%f'))

    ind = np.abs([parse_time(list) - parse_time(img_time) for list in time])
    map_aia = Map(f_aia[ind.argmin(0)])
#    aia_map_updated_pointing = update_pointing(map_aia)
    aia_map_registered = register(map_aia)
#    aia_map_corrected = correct_degradation(aia_map_registered)
    aia_map = aia_map_registered/aia_map_registered.exposure_time
    
# Download the HMI data
    hmi_instr = a.Instrument('HMI')
    hmi_obs = a.Physobs.los_magnetic_field

    f_hmi, file_time = get_hmi_filelist(SDO_data_directory, hmi_obs, file_date)
    if f_hmi == []:
        print('Downloading HMI data')
        data_dir = SDO_data_directory+'HMI_mag/'+file_date
        os.makedirs(data_dir, exist_ok='True')
        result = Fido.search(attrs_time, hmi_instr, hmi_obs, cadence)
        files = Fido.fetch(result, path = data_dir, overwrite=False, progress=True)

    f_hmi, file_time = get_hmi_filelist(SDO_data_directory, hmi_obs, file_date)
    nearest_img = closest(np.array(file_time), dt.datetime.strptime(start_time,"%Y/%m/%d %H:%M:%S"))
    if dt.datetime.strptime(start_time, "%Y/%m/%d %H:%M:%S") < file_time[nearest_img] < dt.datetime.strptime(end_time, "%Y/%m/%d %H:%M:%S"):
        print('HMI Data already downloaded')
    else:
        print('Downloading HMI data')
        data_dir = SDO_data_directory+'HMI_mag/'+file_date
        os.makedirs(data_dir, exist_ok='True')
        result = Fido.search(attrs_time, hmi_instr, hmi_obs, cadence)
        files = Fido.fetch(result, path = data_dir, overwrite=False, progress=True)
        f_hmi = glob.glob(data_dir+'/*.fits')
        f_hmi.sort()

# Now find the nearest HMI image in time to the EIS map and load it as a map
    time = []
    for file_i in f_hmi:
        hdr = fits.getheader(file_i, 1)
        try:
            time.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%fZ'))
        except:
            time.append(dt.datetime.strptime(hdr.get('DATE-OBS'),'%Y-%m-%dT%H:%M:%S.%f'))

    ind = np.abs([parse_time(list) - parse_time(img_time) for list in time])
    map_hmi = Map(f_hmi[ind.argmin(0)])
    hmi_map = map_hmi.reproject_to(aia_map.wcs)

    return aia_map, hmi_map


# %% [markdown]
# Define a function based on Andy's asheis.py class to fit the data for a defined line to get a defined property

# %%
###changed ncpu to 20 hopefully works
def fit_data(file,fitted_lines,line,product,output_location):
    template_name = fitted_lines[f'{line}'][0]
    path = output_location+'/'+file.split('/')[6][0:19]+'.'+(f'{template_name}'.replace(".template",f"-{fitted_lines[f'{line}'][1]}.fit"))
    if os.path.isfile(path) == False:
        template = eispac.read_template(eispac.data.get_fit_template_filepath(template_name))
        cube = eispac.read_cube(file, window=template.central_wave)
        fit_res = eispac.fit_spectra(cube, template, ncpu='30')
        save_filepaths = eispac.save_fit(fit_res, save_dir=output_location)
        print(save_filepaths)
    else:
        fit_res=eispac.read_fit(path)

    fit_res.fit[f'{product}'] = fit_res.shift2wave(fit_res.fit[f'{product}'],wave=195.119)
    return fit_res

# %% [markdown]
# Define a function to get the line intensity

# %%
def get_intensity(line,file,fitted_lines,output_location):
    fit_res = fit_data(file,fitted_lines,line,'int',output_location) # Get fitdata
    m = fit_res.get_map(fitted_lines[f'{line}'][1],measurement='intensity') # From fitdata get map

    save_dir = os.path.join(output_location, 'intensity_files')
    os.makedirs(save_dir, exist_ok=True)
  
    timestamp = m.date.strftime("%Y%m%d_%H%M%S")
    #output_filename = os.path.join(save_dir, f'eis_{m_comp.date.strftime("%Y%m%d_%H%M%S")}_intensity_{linepair}.fits')
    output_filename = os.path.join(save_dir, f'eis_{m.date.strftime("%Y%m%d_%H%M%S")}_intensity_{line}.fits')
    m.save(output_filename, overwrite=True)
    return m, fit_res
    
    
# %% [markdown]
# Define a function to get the Doppler velocity

# %%
def get_velocity(line,file,fitted_lines,output_location):
    fit_res = fit_data(file,fitted_lines,line,'vel',output_location)
    m = fit_res.get_map(component = fitted_lines[f'{line}'][1],measurement='velocity')
    return m, fit_res

# %% [markdown]
# Define a function to get the line width

# %%
def get_width(line,file,fitted_lines,output_location):
    fit_res = fit_data(file,fitted_lines,line,'vel',output_location)
    m = fit_res.get_map(component = fitted_lines[f'{line}'][1],measurement='width')
    return m, fit_res

# %% [markdown]
# Define a function to get and plot the composition using the Si/S line ratio (line ratio here, not full FIP map)

# %%
def get_composition(linepair, filename, output_location, fitted_lines):
    if linepair == "SiS":
        lines=['si_10_258','s_10_264','Si X-S X 1']
    elif linepair == "CaAr":
        lines=['ca_14_193','ar_14_194','Ca XIV-Ar XIV 1'] 
    else:
        print('No line database can be found. Add your line in code.')
        return
        
#    template_names=[fitted_lines[lines[0]][0],fitted_lines[lines[1]][0]]
#    templates = [eispac.read_template(eispac.data.get_fit_template_filepath(t)) for t in template_names]
#
#    outfiles=[]
#    for temp in range(0, len(templates)):
#    
#        t = template_names[temp]
#        path = output_location+'/'+filename.split('/')[5][0:19]+'.'+(f'{t}'.replace(".template",f"-{fitted_lines[lines[temp]][1]}.fit"))
#
#        if os.path.isfile(path) == False:
#            cube = eispac.read_cube(filename, window=t.central_wave)
#            fit_res = eispac.fit_spectra(cube, t, ncpu='max')
#            fit_res.fit['int'] = fit_res.shift2wave(fit_res.fit['int'],wave=195.119)
#        else:
#            fit_res=eispac.read_fit(path)
#    
#        m = fit_res.get_map(component = fitted_lines[lines[temp]][1], measurement='intensity')
#        date = filename.split('/')[5][4:19]
#        os.makedirs(output_location+'/composition_files/', exist_ok=True)
#        m.save(output_location+'/composition_files/eis_'+date+'_'+lines[temp]+'.fits', overwrite=True)
#        outfiles.append(output_location+'/composition_files/eis_'+date+'_'+lines[temp]+'.fits')
#
#    m_upr = Map(outfiles[0])
#    m_lwr = Map(outfiles[1])

    #f_upr = glob.glob(output_location+'/EIS_fit_'+lines[0]+'.asdf')
    #arrs_upr = asdf.open(f_upr[0])
    #m_upr = arrs_upr['int_map']

    f_upr = glob.glob(output_location+'/EIS_fit_'+lines[0]+'.asdf')
    if not f_upr:
        raise FileNotFoundError(f"Upper line file not found: {output_location+'/EIS_fit_'+lines[0]+'.asdf'}")
    arrs_upr = asdf.open(f_upr[0])
    m_upr = arrs_upr['int_map']

    f_lwr = glob.glob(output_location+'/EIS_fit_'+lines[1]+'.asdf')
    if not f_lwr:
        raise FileNotFoundError(f"Lower line file not found: {output_location+'/EIS_fit_'+lines[1]+'.asdf'}")
    arrs_lwr = asdf.open(f_lwr[0])
    m_lwr = arrs_lwr['int_map']
    
    #f_lwr = glob.glob(output_location+'/EIS_fit_'+lines[1]+'.asdf')
    #arrs_lwr = asdf.open(f_lwr[0])
    #m_lwr = arrs_lwr['int_map']

    #m_comp = m_upr
    #m_comp.meta['line_id'] = lines[2]
    #m_comp = Map(m_upr.data/m_lwr.data, m_upr.meta)
    #m_comp.save(output_location+'/composition_files/eis_'+m_comp.date.strftime('%Y%m%d_%H%M%S')+'_composition_'+linepair+'.fits', overwrite=True)
    #return m_comp

    m_comp = Map(m_upr.data/m_lwr.data, m_upr.meta)
    m_comp.meta['line_id'] = lines[2]

    # Ensure the output directory exists before saving
    save_dir = os.path.join(output_location, 'composition_files')
    os.makedirs(save_dir, exist_ok=True)  # Creates the directory if it doesn't exist

    output_filename = os.path.join(save_dir, f'eis_{m_comp.date.strftime("%Y%m%d_%H%M%S")}_composition_{linepair}.fits')

    m_comp.save(output_filename, overwrite=True)
    return m_comp

# %% [markdown]
# Define a function to plot the intensity, Doppler velocity and line width

# %%
def plot_eis_fits(line, int, vel, wid, aia_map, output_location, fitted_lines):
    figs = (12,5)
    wid_rat = [1,1,1,1]
    asp = 1/4
    fig = plt.figure(constrained_layout=True, figsize=figs)
    gs = gridspec.GridSpec(1,4,width_ratios=wid_rat)
    plt.rcParams['font.size'] = '10'
    date = int.date.strftime("%Y%m%d_%H%M%S")
    alpha = 0.1
    plot_name = fitted_lines[f'{line}'][2]+'; '+int.date.strftime("%Y/%m/%dT%H:%M:%S")

# Intensity
    ax1 = fig.add_subplot(gs[0,0], projection = int, label = 'a)')

    lwr_bnd = np.percentile(int.data, alpha)
    upr_bnd = np.percentile(int.data, 100-alpha)
    norm = colors.Normalize(vmin=lwr_bnd, vmax=upr_bnd)
    int.plot_settings['norm'] = norm

    int.plot(axes=ax1, title = 'a) Peak intensity', aspect=asp)
    x = ax1.coords[0]
    x.set_axislabel(' ')
    x.set_ticklabel(exclude_overlapping=True)

    plt.colorbar(ax=ax1,location='right', label='')

# Doppler velocity
    ax2 = fig.add_subplot(gs[0,1], projection = vel, label='b)')

    vel.plot_settings['norm'].vmin = -15
    vel.plot_settings['norm'].vmax = 15
    vel.plot(axes=ax2, title='b) Doppler velocity', aspect=asp)
    x = ax2.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    y = ax2.coords[1]
    y.set_ticklabel_visible(False)

    plt.colorbar(ax=ax2,location='right', label='')

# Line width 
    ax3 = fig.add_subplot(gs[0,2], projection = wid, label='c)')

    lwr_bnd = np.percentile(wid.data, alpha)
    upr_bnd = np.percentile(wid.data, 100-alpha)
    norm = colors.Normalize(vmin=lwr_bnd, vmax=upr_bnd)
    wid.plot_settings['norm'] = norm
    wid.plot(axes=ax3, title='c) Line width', aspect=asp)
    x = ax3.coords[0]
    x.set_axislabel(' ')
    x.set_ticklabel(exclude_overlapping=True)
    y = ax3.coords[1]
    y.set_ticklabel_visible(False)

    plt.colorbar(ax=ax3,location='right', label='')

# AIA context image
    b_left = [int.bottom_left_coord.Tx-200*u.arcsec, int.bottom_left_coord.Ty-200*u.arcsec]
    t_right = [int.top_right_coord.Tx+200*u.arcsec, int.top_right_coord.Ty+200*u.arcsec]

    top_right = SkyCoord(t_right[0],  t_right[1], frame=aia_map.coordinate_frame)
    bottom_left = SkyCoord(b_left[0], b_left[1], frame=aia_map.coordinate_frame)
    aia_map = aia_map.submap(bottom_left, top_right=top_right)

    ax4 = fig.add_subplot(gs[0,3], projection = aia_map, label='d)')

    proc_img = enhance.mgn(aia_map.data, h=0.94, gamma_min=0, gamma_max=aia_map.data.max())
    aia_map = Map(proc_img, aia_map.meta)

    alpha = 0.01
    lwr_bnd = np.percentile(aia_map.data, alpha)
    upr_bnd = np.percentile(aia_map.data, 100-alpha)
    norm = colors.Normalize(vmin=lwr_bnd, vmax=upr_bnd)
    aia_map.plot_settings['norm'] = norm
    aia_map.plot_settings['norm'] = norm
    
    y=ax4.coords[1]
    y.set_axislabel(' ')
    y.set_ticklabel_position('r')
    y.set_axislabel_position('r')

    aia_map.plot(axes=ax4, title='d) AIA context')

#Overplot the EIS FoV
    bottom_left = int.bottom_left_coord
    top_right = int.top_right_coord
    aia_map.draw_quadrangle(bottom_left, top_right=top_right, axes = ax4, edgecolor='blue') 

    plt.suptitle(plot_name)

    plt.savefig(output_location+'/EIS_fits_'+date+'_'+line+'.png', bbox_inches='tight')
    plt.close(fig)

# %% [markdown]
# Define a function to get the composition line ratio

# %%
def plot_intensity(line, intensity_map, aia_map, hmi_map, output_location):
    #date = comp.date.strftime("%Y%m%d_%H%M%S")
    date = intensity_map.date.strftime("%Y%m%d_%H%M%S")

    if intensity_map.dimensions[1]/intensity_map.dimensions[0] >= 20:
        figs = (12,5)
        wid_rat = [1,3,3]
        asp = 1/4
    else:
        figs = (12,6)
        wid_rat = [1,2,2]
        asp = 1/3
    fig = plt.figure(constrained_layout=True, figsize=figs)
    gs = gridspec.GridSpec(1,3,width_ratios=wid_rat)
    plt.rcParams['font.size'] = '10'
    # Intensity
    ax1 = fig.add_subplot(gs[0,0], projection = intensity_map, label = 'a)')
    norm = colors.Normalize(vmin=0, vmax=4)
    comp.plot_settings['norm'] = norm
    comp.plot_settings['cmap'] = 'RdYlBu'
    comp.plot(axes=ax1, title = 'a) '+title, aspect=asp)
    x=ax1.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    plt.colorbar(location='right', label='')
    
    #fig, ax = plt.subplots(figsize=(6,5))
    #intensity_map.plot(axes=ax, title=f'Intensity Map for {line}')
   # plt.colorbar(ax=ax, location='right', label='Intensity')
   # plt.savefig(os.path.join(output_location, f'intensity_map_{line}.png'))

def plot_composition(linepair, comp, aia_map, hmi_map, output_location):
    date = comp.date.strftime("%Y%m%d_%H%M%S")
    if linepair == "SiS":
        title = 'Si/S composition '+date
    elif linepair == "CaAr":
        title = 'Ca/Ar composition '+date

    if comp.dimensions[1]/comp.dimensions[0] >= 20:
        figs = (12,5)
        wid_rat = [1,3,3]
        asp = 1/4
    else:
        figs = (12,6)
        wid_rat = [1,2,2]
        asp = 1/3
    fig = plt.figure(constrained_layout=True, figsize=figs)
    gs = gridspec.GridSpec(1,3,width_ratios=wid_rat)
    plt.rcParams['font.size'] = '10'

# Composition
    ax1 = fig.add_subplot(gs[0,0], projection = comp, label = 'a)')
    norm = colors.Normalize(vmin=0, vmax=4)
    comp.plot_settings['norm'] = norm
    comp.plot_settings['cmap'] = 'RdYlBu'
    comp.plot(axes=ax1, title = 'a) '+title, aspect=asp)
    x=ax1.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    plt.colorbar(location='right', label='')
    
# Process the HMI and AIA context images
    b_left = [comp.bottom_left_coord.Tx-200*u.arcsec, comp.bottom_left_coord.Ty-200*u.arcsec]
    t_right = [comp.top_right_coord.Tx+200*u.arcsec, comp.top_right_coord.Ty+200*u.arcsec]

    top_right = SkyCoord(t_right[0],  t_right[1], frame=aia_map.coordinate_frame)
    bottom_left = SkyCoord(b_left[0], b_left[1], frame=aia_map.coordinate_frame)
    aia_map = aia_map.submap(bottom_left, top_right=top_right)
    hmi_map = hmi_map.submap(bottom_left, top_right=top_right)

# AIA context image
    ax2 = fig.add_subplot(gs[0,1], projection = aia_map, label='b)')

    proc_img = enhance.mgn(aia_map.data, h=0.94, gamma_min=0, gamma_max=aia_map.data.max())
    aia_map = Map(proc_img, aia_map.meta)

    alpha = 0.01
    lwr_bnd = np.percentile(aia_map.data, alpha)
    upr_bnd = np.percentile(aia_map.data, 100-alpha)
    norm = colors.Normalize(vmin=lwr_bnd, vmax=upr_bnd)
    aia_map.plot_settings['norm'] = norm
    aia_map.plot_settings['norm'] = norm
    
    aia_map.plot(axes=ax2, title='b) AIA context')
    y = ax2.coords[1]
    y.set_axislabel(' ')
    y.set_ticklabel_position('l')
    y.set_axislabel_position('l')

#Overplot the EIS FoV
    bottom_left = comp.bottom_left_coord
    top_right = comp.top_right_coord
    aia_map.draw_quadrangle(bottom_left, top_right=top_right, axes = ax2, edgecolor='blue') 

# HMI context image
    ax3 = fig.add_subplot(gs[0,2], projection = hmi_map, label='c)')

    norm = colors.Normalize(vmin=-200, vmax=200)
    hmi_map.plot_settings['norm'] = norm
    
    hmi_map.plot(axes=ax3, title='c) HMI context')
    y=ax3.coords[1]
    y.set_axislabel(' ')
    y.set_ticklabel_visible(False)


#Overplot the EIS FoV
    bottom_left = comp.bottom_left_coord
    top_right = comp.top_right_coord
    hmi_map.draw_quadrangle(bottom_left, top_right=top_right, axes = ax3, edgecolor='blue') 

    plt.savefig(output_location+'/EIS_composition_'+date+'_'+linepair+'.png', bbox_inches='tight')
    plt.close(fig)

# %% [markdown]
# Use the previously defined functions to process and plot the data

def run_eis_processing():

    eis_evts = ['20151018_173743']
    time_range = ['2015/10/18T17:36:00', '2015/10/18T17:38:00']
    date = '20151018'

    timerange = [dt.datetime.strptime(tr,'%Y/%m/%dT%H:%M:%S') for tr in time_range]
    file_date = dt.datetime.strftime(dt.datetime.strptime(time_range[0],'%Y/%m/%dT%H:%M:%S'), '%Y/%m/%d')

    data_location = '/mnt/scratch/data/spruksk2/data'
    os.makedirs(data_location, exist_ok=True)
    output_location = '/mnt/scratch/data/spruksk2/python_output/EIS_work/'+date
    os.makedirs(output_location, exist_ok=True)

    f_eis = glob.glob(data_location+'/*.data.h5')
    f_eis.sort()

    #results = Fido.search(a.Time(timerange[0], timerange[1]),
    #                      a.Instrument('EIS'),
    #                      a.Physobs.intensity,
    #                      a.Source('Hinode'),
    #                      a.Provider('NRL'),
    #                      a.Level('1'))  
    #files = Fido.fetch(results, path = data_location, overwrite=False, progress=True)

# %%
    fitted_lines = {
        "fe_12_195" : ["fe_12_195_119.2c.template.h5",0,"Fe XII 195$\AA$"],
        "ar_14_194" : ["ar_14_194_396.6c.template.h5",5,"Ar XIV 194$\AA$"],
        "ca_14_193" : ["ca_14_193_874.6c.template.h5",1,"Ca XIV 193$\AA$"],
        "si_10_258" : ["si_10_258_375.1c.template.h5",0,"Si X 258$\AA$"],
        "s_10_264" : ["s__10_264_233.1c.template.h5",0,"S X 264$\AA$"],
        "fe_13_202" : ["fe_13_202_044.1c.template.h5",0,"Fe XIII 202$\AA$"],
        "fe_13_203" : ["fe_13_203_826.2c.template.h5",1,"Fe XIII 203$\AA$"]
    }


# %%
    #line_list = ["ar_14_194","ca_14_193"]#,"si_10_258","s_10_264"]#,"fe_12_195","fe_13_202","fe_13_203"]
    line_list = ["fe_12_195", "ar_14_194","ca_14_193","si_10_258"]#"s_10_264"]#,"fe_12_195","fe_13_202","fe_13_203"]
    #line_list = ["ca_14_193"]#,"si_10_258","s_10_264"]#,"fe_12_195","fe_13_202","fe_13_203"]

    for event in range(0, len(eis_evts)):

        img_time = dt.datetime.strftime(dt.datetime.strptime(eis_evts[event],'%Y%m%d_%H%M%S'), '%Y/%m/%dT%H:%M:%S')
        date = dt.datetime.strftime(dt.datetime.strptime(eis_evts[event],'%Y%m%d_%H%M%S'), '%Y%m%d')
        file_date = dt.datetime.strftime(dt.datetime.strptime(eis_evts[event],'%Y%m%d_%H%M%S'), '%Y/%m/%d')

        data_location = '/mnt/scratch/data/spruksk2/data'
        os.makedirs(data_location, exist_ok=True)
        output_location = '/mnt/scratch/data/spruksk2/python_output/EIS_work/'+eis_evts[event]
        os.makedirs(output_location+'/save_files', exist_ok=True)
        os.makedirs(output_location+'/plots', exist_ok=True)

        f_eis = glob.glob(data_location+'/eis_'+eis_evts[event]+'.data.h5')
        f_eis.sort()
    
        for file in f_eis:

            print(' ')
            print(file)
            print(' ')

            aia_map, hmi_map = get_context_img(img_time)

            for wvl in line_list:

                print(' ')
                print(wvl)
                print(' ')

                i_map, fit_res = get_intensity(wvl,file,fitted_lines,output_location+'/save_files')
                v_map, fit_res = get_velocity(wvl,file,fitted_lines,output_location+'/save_files')
                w_map, fit_res = get_width(wvl,file,fitted_lines,output_location+'/save_files')
                plot_eis_fits(wvl, i_map, v_map, w_map, aia_map, output_location+'/plots', fitted_lines)

                tree = {'int_map':i_map, 'dopp_map':v_map, 'wid_map':w_map}
                with asdf.AsdfFile(tree) as asdf_file:  
                    asdf_file.write_to(output_location+'/EIS_fit_'+wvl+'.asdf', all_array_compression='zlib')

         
#           m_SiS = get_composition('SiS', file, output_location+'/save_files', fitted_lines)
#           plot_composition('SiS', m_SiS, aia_map, hmi_map, output_location+'/plots')
            #m_CaAr = get_composition('CaAr', file, output_location, fitted_lines)
            #plot_composition('CaAr', m_CaAr, aia_map, hmi_map, output_locationputty+'/plots')

            # Save intensity map in FITS format (similar to CaAr)
            #intensity_filename = os.path.join(output_location, 'intensity_files', f'eis_{wvl}_intensity.fits')
            #i_map.save(intensity_filename, overwrite=True)

            # Optional: Convert intensity map to SunPy map and save again
            m_intensity = i_map
            plot_intensity(wvl, m_intensity, aia_map, hmi_map, output_location+'/plots')
            #plot_eis_fits(wvl, m_intensity, aia_map, hmi_map, output_location+ '/plots') 
            

if __name__ == "__main__":
    run_eis_processing()
