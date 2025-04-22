import asdf
from sunpy.map import Map
import glob
import matplotlib.pyplot as mpl
import numpy as np
from reproject.mosaicking import reproject_and_coadd, find_optimal_celestial_wcs
from reproject import reproject_adaptive

f_0 = glob.glob('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/EIS_fit_fe_12_195.asdf')
arrs = asdf.open(f_0[0])
dopp_map_0 = arrs['dopp_map']

f_1 = glob.glob('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/EIS_fit_fe_12_195.asdf')
arrs = asdf.open(f_1[0])
dopp_map_1 = arrs['dopp_map']

f_2 = glob.glob('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/EIS_fit_fe_12_195.asdf')
arrs = asdf.open(f_2[0])
dopp_map_2 = arrs['dopp_map']

f_3 = glob.glob('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/EIS_fit_fe_12_195.asdf')
arrs = asdf.open(f_3[0])
dopp_map_3 = arrs['dopp_map']

wcs_out, shape_out = find_optimal_celestial_wcs(comp_array)

array, footprint = reproject_and_coadd(comp_array, wcs_out, shape_out, reproject_function=reproject_adaptive)

combined_map = Map(array, wcs_out)

aspect_ratio = combined_map.meta['CDELT2'] / np.abs(combined_map.meta['CDELT1'])
combined_map.peek(vmin=-20, vmax=20, cmap=mpl.colormaps['coolwarm'], aspect=aspect_ratio)


