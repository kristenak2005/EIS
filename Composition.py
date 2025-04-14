import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import matplotlib.colora as colors


i_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/composition_files/eis_20151018_102719_composition_CaAr.fits")   
i_map1 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/composition_files/eis_20151018_113839_composition_CaAr.fits") 
i_map2 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/composition_files/eis_20151018_124939_composition_CaAr.fits")
i_map3 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_140939/composition_files/eis_20151018_140939_composition_CaAr.fits")
i_map4 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_160113/composition_files/eis_20151018_160113_composition_CaAr.fits")
i_map5 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_173743/composition_files/eis_20151018_173743_composition_CaAr.fits")
i_map6 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/composition_files/eis_20151018_191443_composition_CaAr.fits")
i_map7 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/composition_files/eis_20151018_205143_composition_CaAr.fits")
i_map8 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_222743/composition_files/eis_20151018_222743_composition_CaAr.fits")
i_map9 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_235413/composition_files/eis_20151018_235413_composition_CaAr.fits")
i_map10 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_012613/composition_files/eis_20151019_012613_composition_CaAr.fits")
i_map11 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_030043/composition_files/eis_20151019_030043_composition_CaAr.fits")
i_map12 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_043942/composition_files/eis_20151019_043942_composition_CaAr.fits")
i_map13 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_061712/composition_files/eis_20151019_061712_composition_CaAr.fits")
i_map14 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_075442/composition_files/eis_20151019_075442_composition_CaAr.fits")
i_map15 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_093012/composition_files/eis_20151019_093012_composition_CaAr.fits")
i_map16 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_105442/composition_files/eis_20151019_105442_composition_CaAr.fits")
i_map17 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_120542/composition_files/eis_20151019_120542_composition_CaAr.fits")
i_map18 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_131642/composition_files/eis_20151019_131642_composition_CaAr.fits")
i_map19 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_144842/composition_files/eis_20151019_144842_composition_CaAr.fits")
i_map20 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_163612/composition_files/eis_20151019_163612_composition_CaAr.fits")

#v_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_velocity_fe_12_195.fits")
#w_map = Map('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_width_fe_12_195.fits')

bottom_left = [200, -800] * u.arcsec
top_right = [1000, -200] * u.arcsec

#i_smap = i_map.submap(SkyCoord(*bottom_left, frame=i_map.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map.coordinate_frame))

#i_smap1 = i_map1.submap(SkyCoord(*bottom_left, frame=i_map1.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map1.coordinate_frame))

#i_smap2 = i_map2.submap(SkyCoord(*bottom_left, frame=i_map2.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map2.coordinate_frame))

#i_smap3 = i_map3.submap(SkyCoord(*bottom_left, frame=i_map3.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map3.coordinate_frame))

#i_smap4 = i_map4.submap(SkyCoord(*bottom_left, frame=i_map4.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map4.coordinate_frame))

#i_smap5 = i_map5.submap(SkyCoord(*bottom_left, frame=i_map5.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map5.coordinate_frame))

#i_smap6 = i_map6.submap(SkyCoord(*bottom_left, frame=i_map6.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map6.coordinate_frame))

#i_smap7 = i_map7.submap(SkyCoord(*bottom_left, frame=i_ma7p.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map7.coordinate_frame))

#i_smap8 = i_map8.submap(SkyCoord(*bottom_left, frame=i_map8.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map8.coordinate_frame))

#i_smap9 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map9.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map9.coordinate_frame))

#i_smap10 = i_map10.submap(SkyCoord(*bottom_left, frame=i_map10.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map10.coordinate_frame))

#i_smap11 = i_map11.submap(SkyCoord(*bottom_left, frame=i_map11.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map11.coordinate_frame))

#i_smap12 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map12.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map12.coordinate_frame))

#i_smap13 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map13.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map13.coordinate_frame)) 

#i_smap14 = i_map14.submap(SkyCoord(*bottom_left, frame=i_map14.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map14.coordinate_frame))

#i_smap15 = i_map15.submap(SkyCoord(*bottom_left, frame=i_map15.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map15.coordinate_frame))

#i_smap16 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map16.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map16.coordinate_frame))

#i_smap17 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map17.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map17.coordinate_frame))

#i_smap18 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map18.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map18.coordinate_frame))

#i_smap19 = i_map19.submap(SkyCoord(*bottom_left, frame=i_map19.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map19.coordinate_frame))

#i_smap20 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map20.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map20.coordinate_frame))




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
i_map10.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map11.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map12.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map13.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map14.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map15.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map16.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map17.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map18.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map19.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')
i_map20.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='RdYlBu')

aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Composition Composite Map")
plt.show()
