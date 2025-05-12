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
output_location_ar = '/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_205143/'
file_ar = glob.glob(output_location_ar + '/eis_20151018_205143_composition_FeS.fits')
compFeS = Map(file_ar[0])  # Existing CaAr composition map

ar_bottom_left = SkyCoord(300 * u.arcsec, -550 * u.arcsec, frame=compFeS.coordinate_frame)
ar_top_right = SkyCoord(700 * u.arcsec, -250 * u.arcsec, frame=compFeS.coordinate_frame)
comp_ar_map = compFeS.submap(ar_bottom_left, top_right=ar_top_right)
comp_ar = comp_ar_map.data.flatten()

#Quiet sun
qs_bottom_left = SkyCoord(100 * u.arcsec, 0 * u.arcsec, frame=compFeS.coordinate_frame)
qs_top_right = SkyCoord(500 * u.arcsec, 300 * u.arcsec, frame=compFeS.coordinate_frame)
comp_qs_map = compCaAr.submap(qs_bottom_left, top_right=qs_top_right)
comp_qs = comp_qs_map.data.flatten()

output_location_qs = '/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_191443'
file_qs = glob.glob(output_location_qs + '/eis_20151018_191443_composition_FeS.fits')
compFeS_qs = Map(file_qs[0])

qs_map_file = compFeS_qs.submap(qs_bottom_left, top_right=qs_top_right)
comp_qs_file = qs_map_file.data.flatten()

ar_data = pd.DataFrame({'comp': comp_ar, 'type': np.tile('Active Region', len(comp_ar))})
qs_data = pd.DataFrame({'comp': comp_qs, 'type': np.tile('Quiet Sun (same map)', len(comp_qs))})
qs_file_data = pd.DataFrame({'comp': comp_qs_file, 'type': np.tile('Quiet Sun (separate file)', len(comp_qs_file))})

# Combine all
df = pd.concat([ar_data, qs_data, qs_file_data], ignore_index=True)

# --- Plot KDE ---
fig = plt.figure(figsize=(15,14))
sns.kdeplot(df, x='comp', hue='type', fill=True, alpha=0.5, legend=True, bw_adjust=0.5)
plt.show()
