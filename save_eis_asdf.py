import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u
import asdf

eis_map =sunpy.map.Map("eis_20151018_102719_sis.fits")

eis_map.peek()
