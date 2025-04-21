import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord


i_map = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_102719/eis_20151018_102719_FeS.fits")
i_map1 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_160113/eis_20151018_160113_FeS.fits")
i_map2 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_173743/eis_20151018_173743_FeS.fits")
i_map3 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_191443/eis_20151018_191443_FeS.fits")
i_map4 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_205143/eis_20151018_205143_FeS.fits")
i_map5 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151018_222743/eis_20151018_222743_FeS.fits")
i_map6 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_012613/eis_20151019_1012613_FeS.fits")
i_map7 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_030043/eis_20151019_030043_FeS.fits")
i_map8 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_043942/eis_20151019_043942_FeS.fits")
i_map9 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_061712/eis_20151019_161712_FeS.fits")
i_map10 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_075442/eis_20151019_075442_FeS.fits")
i_map11 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_120542/eis_20151019_120542_FeS.fits")
i_map12 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_131642/eis_20151019_131642_FeS.fits")
i_map13 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_144842/eis_20151019_144842_FeS.fits")
i_map14 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_163612/eis_20151019_163612_FeS.fits")
i_map15 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_212642/eis_20151019_212642_FeS.fits")
i_map16 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151019_230142/eis_20151019_230142_FeS.fits")
i_map17 = Map("/mnt/scratch/data/spruksk2/extracted_files/eis_20151020_002612/eis_20151020_002612_FeS.fits")

bottom_left = [200, -800] * u.arcsec
top_right = [1000, -200] * u.arcsec

fig = plt.figure()
ax = fig.add_subplot(projection=i_map)

i_map.plot(axes=ax, clip_interval=(1, 99.9) * u.percent, zorder=0,autoalign=True, cmap='RdYlBu')
i_map1.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=1, autoalign=True, cmap='RdYlBu')
i_map2.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=2, autoalign=True, cmap='RdYlBu')
i_map3.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=3, autoalign=True, cmap='RdYlBu')
i_map4.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=4, autoalign=True, cmap='RdYlBu')
i_map5.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=5, autoalign=True, cmap='RdYlBu')
i_map6.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=6, autoalign=True, cmap='RdYlBu')
i_map7.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=7, autoalign=True, cmap='RdYlBu')
i_map8.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=8, autoalign=True, cmap='RdYlBu')
i_map9.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map10.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=10, autoalign=True, cmap='RdYlBu')
i_map11.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=11, autoalign=True, cmap='RdYlBu')
i_map12.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=12, autoalign=True, cmap='RdYlBu')
i_map13.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=13, autoalign=True, cmap='RdYlBu')
i_map14.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=14, autoalign=True, cmap='RdYlBu')
i_map15.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=15, autoalign=True, cmap='RdYlBu')
i_map16.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=16, autoalign=True, cmap='RdYlBu')
i_map17.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=17, autoalign=True, cmap='RdYlBu')

aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Composition Composite Map")
plt.show()
