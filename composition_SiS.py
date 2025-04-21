import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord


i_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/composition_files/eis_20151018_102719_composition_SiS.fits")   
i_map1 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/composition_files/eis_20151018_113839_composition_SiSfits") 
i_map2 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/composition_files/eis_20151018_124939_composition_SiS.fits")
i_map3 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_140939/composition_files/eis_20151018_140939_composition_SiS.fits")
i_map4 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_160113/composition_files/eis_20151018_160113_composition_SiS.fits")
i_map5 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_173743/composition_files/eis_20151018_173743_composition_SiS.fits")
i_map6 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/composition_files/eis_20151018_191443_composition_SiS.fits")
i_map7 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/composition_files/eis_20151018_205143_composition_SiS.fits")
i_map8 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_222743/composition_files/eis_20151018_222743_composition_SiS.fits")
i_map9 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_235413/composition_files/eis_20151018_235413_composition_SiS.fits")
i_map10 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_012613/composition_files/eis_20151019_012613_composition_SiS.fits")
i_map11 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_030043/composition_files/eis_20151019_030043_composition_SiS.fits")
i_map12 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_043942/composition_files/eis_20151019_043942_composition_SiS.fits")
i_map13 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_061712/composition_files/eis_20151019_061712_composition_SiS.fits")
i_map14 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_075442/composition_files/eis_20151019_075442_composition_SiS.fits")
i_map15 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_093012/composition_files/eis_20151019_093012_composition_SiS.fits")
i_map16 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_105442/composition_files/eis_20151019_105442_composition_SiS.fits")
i_map17 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_120542/composition_files/eis_20151019_120542_composition_SiS.fits")
i_map18 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_131642/composition_files/eis_20151019_131642_composition_SiS.fits")
i_map19 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_144842/composition_files/eis_20151019_144842_composition_SiS.fits")
i_map20 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_163612/composition_files/eis_20151019_163612_composition_SiS.fits")
i_map21 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_181312/composition_files/eis_20151019_181312_composition_SiS.fits")
i_map22 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_194942/composition_files/eis_20151019_194942_composition_SiS.fits")
i_map23 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_212642/composition_files/eis_20151019_212642_composition_SiS.fits")
i_map24 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_230142/composition_files/eis_20151019_230142_composition_SiS.fits")
i_map25 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151020_002612/composition_files/eis_20151020_002612_composition_SiS.fits")


#v_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_velocity_fe_12_195.fits")
#w_map = Map('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_width_fe_12_195.fits')

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
i_map18.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=18, autoalign=True, cmap='RdYlBu')
i_map19.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=19, autoalign=True, cmap='RdYlBu')
i_map20.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=20, autoalign=True, cmap='RdYlBu')
i_map21.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=21, autoalign=True, cmap='RdYlBu')
i_map22.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=22, autoalign=True, cmap='RdYlBu')
i_map23.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=23, autoalign=True, cmap='RdYlBu')
i_map24.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=24, autoalign=True, cmap='RdYlBu')
i_map25.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=25, autoalign=True, cmap='RdYlBu')


aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Composition Composite Map")
plt.show()
