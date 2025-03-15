from sunpy.map import Map
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.coordinates import SphericalScreen

# Load your EIS FITS maps
i_map = Map('eis_intensity.fits')  #Example: eis_20151018_102719_intensity_fe_12_195.fits
v_map = Map('eis_velocity.fits') #Example: eis_20151018_102719_velocity_fe_12_195.fits

w_map = Map('eis_width.fits')


with SphericalScreen(i_map.observer_coordinate):
    v_map = v_map.reproject_to(i_map.wcs)
    w_map = w_map.reproject_to(i_map.wcs)
