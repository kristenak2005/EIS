import glob
import asdf
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sunpy.coordinates import propagate_with_solar_surface

output_location = '/cd/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/composition_files'
#output_location1 = '/cd/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/composition_files'

file = glob.glob(output_location+'/eis_20151018_205143_composition_CaAr.fits')

arrs = asdf.open(file[0])

comp_map = arrs['int_map']
#int_map = arrs['int_map']
#dopp_map = arrs['dopp_map']
#wid_map = arrs['wid_map']

#active sun
ar_bottom_left = SkyCoord(300 * u.arcsec, -550 * u.arcsec, frame=comp_map.coordinate_frame)
ar_top_right = SkyCoord(700 * u.arcsec, -250 * u.arcsec, frame=comp_map.coordinate_frame)

comp_ar_map = comp_map.submap(ar_bottom_left, top_right = ar_top_right)
#int_ar_map = int_map.submap(ar_bottom_left, top_right = ar_top_right)
#dopp_ar_map = dopp_map.submap(ar_bottom_left, top_right = ar_top_right)
#wid_ar_map = wid_map.submap(ar_bottom_left, top_right = ar_top_right)

comp_ar = comp_ar_map.data.flatten()
#int_ar = int_ar_map.data.flatten()
#dopp_ar = dopp_ar_map.data.flatten()
#wid_ar = wid_ar_map.data.flatten()

#quiet sun
qs_bottom_left = SkyCoord(300 * u.arcsec, -550 * u.arcsec, frame=comp_map.coordinate_frame)
qs_top_right = SkyCoord(650 * u.arcsec, -250 * u.arcsec, frame=comp_map.coordinate_frame)

comp_qs_map = comp_map.submap(qs_bottom_left, top_right = qs_top_right)
#int_qs_map = int_map.submap(qs_bottom_left, top_right = qs_top_right)
#dopp_qs_map = dopp_map.submap(qs_bottom_left, top_right = qs_top_right)
#wid_qs_map = wid_map.submap(qs_bottom_left, top_right = qs_top_right)

comp_qs = comp_qs_map.data.flatten()
#int_qs = int_qs_map.data.flatten()
#dopp_qs = dopp_qs_map.data.flatten()
#wid_qs = wid_qs_map.data.flatten()


#ar_bottom_left = SkyCoord(500 * u.arcsec, -60 * u.arcsec, frame=int_map.coordinate_frame)
#ar_top_right = SkyCoord(560 * u.arcsec, 150 * u.arcsec, frame=int_map.coordinate_frame)

#with propagate_with_solar_surface():
        #diffrot_ar_bl = ar_bottom_left.transform_to(int_map_0.coordinate_frame)

#if (diffrot_ar_bl.Tx.arcsecond < int_map_0.meta['crval1']):
#        diffrot_ar_bl.Tx.arcsecond = int_map_0.meta['crval1']
        

ar_data = pd.DataFrame({'comp':comp_ar, 'type':np.tile('Active Region',len(int_ar))})
#ar_data = pd.DataFrame({'int':int_ar, 'vdopp':dopp_ar, 'wid':wid_ar, 'type':np.tile('Active Region',len(int_ar))})
qs_data = pd.DataFrame({'comp':comp_qs, 'type':np.tile('Quiet Sun',len(int_qs))})

#qs_data = pd.DataFrame({'int':int_qs, 'vdopp':dopp_qs, 'wid':wid_qs, 'type':np.tile('Quiet Sun',len(int_qs))})
df = pd.concat([ar_data, qs_data], ignore_index=True)

fig = plt.figure(figsize=(15,14))

h0=sns.kdeplot(df, x='wid', hue='type', fill = True, alpha = 0.5, legend = True, bw_adjust=.5)

plt.show() 
