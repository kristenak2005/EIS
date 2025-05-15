import glob
import asdf
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sunpy.map import Map

#Active region
output_location_ar = '/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/composition_files'
file_ar = glob.glob(output_location_ar + '/eis_20151018_205143_composition_CaAr.fits')
compCaAr = Map(file_ar[0])  # Existing CaAr composition map

ar_bottom_left = SkyCoord(300 * u.arcsec, -550 * u.arcsec, frame=compCaAr.coordinate_frame)
ar_top_right = SkyCoord(700 * u.arcsec, -250 * u.arcsec, frame=compCaAr.coordinate_frame)
comp_ar_map = compCaAr.submap(ar_bottom_left, top_right=ar_top_right)
comp_ar = comp_ar_map.data.flatten()

#Quiet sun
qs_bottom_left = SkyCoord(100 * u.arcsec, 0 * u.arcsec, frame=compCaAr.coordinate_frame)
qs_top_right = SkyCoord(500 * u.arcsec, 300 * u.arcsec, frame=compCaAr.coordinate_frame)
comp_qs_map = compCaAr.submap(qs_bottom_left, top_right=qs_top_right)
comp_qs = comp_qs_map.data.flatten()

output_location_qs = '/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/composition_files'
file_qs = glob.glob(output_location_qs + '/eis_20151018_191443_composition_CaAr.fits')
compCaAr_qs = Map(file_qs[0])

qs_map_file = compCaAr_qs.submap(qs_bottom_left, top_right=qs_top_right)
comp_qs_file = qs_map_file.data.flatten()

ar_data = pd.DataFrame({'comp': comp_ar, 'type': np.tile('Active Region', len(comp_ar))})
qs_data = pd.DataFrame({'comp': comp_qs, 'type': np.tile('Quiet Sun (same map)', len(comp_qs))})
qs_file_data = pd.DataFrame({'comp': comp_qs_file, 'type': np.tile('Quiet Sun (separate file)', len(comp_qs_file))})

# Combine all
df = pd.concat([ar_data, qs_data, qs_file_data], ignore_index=True)
df.loc[(df['comp'] <= 0) | (df['comp'] > 10), 'comp'] = np.nan

# --- Plot KDE ---
fig = plt.figure(figsize=(15,14))
sns.kdeplot(df, x='comp', hue='type', fill=True, alpha=0.5, legend=True, bw_adjust=0.5)
#plt.xlim(-2000, 2000)
#data[data <= 0] = 0
plt.show()

