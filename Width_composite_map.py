import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord


i_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/width_files/eis_20151018_102719_width_fe_12_195.fits")   
i_map1 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/save_files/width_files/eis_20151018_113839_width_fe_12_195.fits") 
i_map2 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/save_files/width_files/eis_20151018_124939_width_fe_12_195.fits")
i_map3 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_140939/save_files/width_files/eis_20151018_140939_width_fe_12_195.fits")
i_map4 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_160113/save_files/width_files/eis_20151018_160113_width_fe_12_195.fits")
i_map5 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_173743/save_files/width_files/eis_20151018_173743_width_fe_12_195.fits")
i_map6 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/save_files/width_files/eis_20151018_191443_width_fe_12_195.fits")
i_map7 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/save_files/width_files/eis_20151018_205143_width_fe_12_195.fits")
i_map8 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_222743/save_files/width_files/eis_20151018_222743_width_fe_12_195.fits")
i_map9 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_235413/save_files/width_files/eis_20151018_235413_width_fe_12_195.fits")
i_map10 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_012613/save_files/width_files/eis_20151019_012613_width_fe_12_195.fits")
i_map11 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_030043/save_files/width_files/eis_20151019_030043_width_fe_12_195.fits")
i_map12 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_043942/save_files/width_files/eis_20151019_043942_width_fe_12_195.fits")
i_map13 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_061712/save_files/width_files/eis_20151019_061712_width_fe_12_195.fits")
i_map14 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_075442/save_files/width_files/eis_20151019_075442_width_fe_12_195.fits")
i_map15 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_093012/save_files/width_files/eis_20151019_093012_width_fe_12_195.fits")
i_map16 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_105442/save_files/width_files/eis_20151019_105442_width_fe_12_195.fits")
i_map17 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_120542/save_files/width_files/eis_20151019_120542_width_fe_12_195.fits")

#v_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_velocity_fe_12_195.fits")
#w_map = Map('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_width_fe_12_195.fits')

bottom_left = [200, -800] * u.arcsec
top_right = [1000, -200] * u.arcsec

#i_smap = i_map.submap(SkyCoord(*bottom_left, frame=i_map.coordinate_frame),
                  #    top_right=SkyCoord(*top_right, frame=i_map.coordinate_frame))

#i_smap1 = i_map1.submap(SkyCoord(*bottom_left, frame=i_map1.coordinate_frame),
                   #   top_right=SkyCoord(*top_right, frame=i_map1.coordinate_frame))

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


##i_smap11 = i_map11.submap(SkyCoord(*bottom_left, frame=i_map11.coordinate_frame),
  #                    top_right=SkyCoord(*top_right, frame=i_map11.coordinate_frame))


#i_smap12 = i_map12.submap(SkyCoord(*bottom_left, frame=i_map12.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map12.coordinate_frame))


#i_smap13 = i_map13.submap(SkyCoord(*bottom_left, frame=i_map13.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map13.coordinate_frame))


#i_smap14 = i_map14.submap(SkyCoord(*bottom_left, frame=i_map14.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map14.coordinate_frame))


#i_smap15 = i_map15.submap(SkyCoord(*bottom_left, frame=i_map15.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map15.coordinate_frame))


#i_smap16 = i_map16.submap(SkyCoord(*bottom_left, frame=i_map16.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map16.coordinate_frame))


#i_smap17 = i_map17.submap(SkyCoord(*bottom_left, frame=i_map17.coordinate_frame),
 #                     top_right=SkyCoord(*top_right, frame=i_map17.coordinate_frame))



min_value = 0
max_value = 0.4

fig = plt.figure()
ax = fig.add_subplot(projection=i_map)
i_map.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, zorder=0,autoalign=True, cmap='coolwarm')
i_map1.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=1, autoalign=True, cmap='coolwarm')
i_map2.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=2, autoalign=True, cmap='coolwarm')
i_map3.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=3, autoalign=True, cmap='coolwarm')
i_map4.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=4, autoalign=True, cmap='coolwarm')
i_map5.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=5, autoalign=True, cmap='coolwarm')
i_map6.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=6, autoalign=True, cmap='coolwarm')
i_map7.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=7, autoalign=True, cmap='coolwarm')
i_map8.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=8, autoalign=True, cmap='coolwarm')
i_map9.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=9, autoalign=True, cmap='coolwarm')
i_map10.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=10, autoalign=True, cmap='coolwarm')
i_map11.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=11, autoalign=True, cmap='coolwarm')
i_map12.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=12, autoalign=True, cmap='coolwarm')
i_map13.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=13, autoalign=True, cmap='coolwarm')
i_map14.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=14, autoalign=True, cmap='coolwarm')
i_map15.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=15, autoalign=True, cmap='coolwarm')
i_map16.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=16, autoalign=True, cmap='coolwarm')
i_map17.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=17, autoalign=True, cmap='coolwarm')
i_map18.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=18, autoalign=True, cmap='coolwarm')
i_map19.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=19, autoalign=True, cmap='coolwarm')
i_map20.plot(axes=ax, clip_interval=(min_value, max_value) * u.percent, alpha=0.7, zorder=20, autoalign=True, cmap='coolwarm')

aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Width Composite Map")
plt.show()
