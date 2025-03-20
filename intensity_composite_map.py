import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord


i_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/intensity_files/eis_20151018_102719_intensity_fe_12_195.fits")   
i_map1 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/save_files/intensity_files/eis_20151018_113839_intensity_fe_12_195.fits") 
i_map2 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/save_files/intensity_files/eis_20151018_124939_intensity_fe_12_195.fits")
i_map3 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_140939/save_files/intensity_files/eis_20151018_140939_intensity_fe_12_195.fits")
i_map4 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_160113/save_files/intensity_files/eis_20151018_160113_intensity_fe_12_195.fits")
i_map5 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_173743/save_files/intensity_files/eis_20151018_173743_intensity_fe_12_195.fits")
i_map6 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/save_files/intensity_files/eis_20151018_191443_intensity_fe_12_195.fits")
i_map7 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/save_files/intensity_files/eis_20151018_205143_intensity_fe_12_195.fits")
i_map8 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_222743/save_files/intensity_files/eis_20151018_222743_intensity_fe_12_195.fits")
i_map9 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_235413/save_files/intensity_files/eis_20151018_235413_intensity_fe_12_195.fits")




#v_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_velocity_fe_12_195.fits")
#w_map = Map('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_width_fe_12_195.fits')

bottom_left = [200, -800] * u.arcsec
top_right = [1000, -200] * u.arcsec

i_smap = i_map.submap(SkyCoord(*bottom_left, frame=i_map.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map.coordinate_frame))

i_smap1 = i_map1.submap(SkyCoord(*bottom_left, frame=i_map1.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map1.coordinate_frame))

i_smap2 = i_map2.submap(SkyCoord(*bottom_left, frame=i_map2.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map2.coordinate_frame))

i_smap3 = i_map3.submap(SkyCoord(*bottom_left, frame=i_map3.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map3.coordinate_frame))

i_smap4 = i_map4.submap(SkyCoord(*bottom_left, frame=i_map4.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map4.coordinate_frame))

i_smap5 = i_map5.submap(SkyCoord(*bottom_left, frame=i_map5.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map5.coordinate_frame))

i_smap6 = i_map6.submap(SkyCoord(*bottom_left, frame=i_map6.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map6.coordinate_frame))

i_smap7 = i_map7.submap(SkyCoord(*bottom_left, frame=i_ma7p.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map7.coordinate_frame))

i_smap8 = i_map8.submap(SkyCoord(*bottom_left, frame=i_map8.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map8.coordinate_frame))

i_smap9 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map9.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map9.coordinate_frame))




fig = plt.figure()
ax = fig.add_subplot(projection=i_map)
i_map.plot(axes=ax, clip_interval=(1, 99.9) * u.percent, zorder=0)
i_map1.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=1)
i_map2.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=2)
i_map3.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=3)
i_map4.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=4)
i_map5.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=5)
i_map6.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=6)
i_map7.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=7)
i_map8.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=8)
i_map9.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=9)
ax.set_title("EIS Intensity Composite Map")
plt.show()
