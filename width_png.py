import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord


i_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/plots/EIS_Width_20151018_102719_fe_12_195.png")   
i_map1 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/plots/EIS_Width_20151018_113839_fe_12_195.png") 
i_map2 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/plots/Eis_Width_20151018_124939_fe_12_195.png")

bottom_left = [200, -800] * u.arcsec
top_right = [1000, -200] * u.arcsec

i_smap = i_map.submap(SkyCoord(*bottom_left, frame=i_map.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map.coordinate_frame))

i_smap1 = i_map1.submap(SkyCoord(*bottom_left, frame=i_map1.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map1.coordinate_frame))

i_smap2 = i_map2.submap(SkyCoord(*bottom_left, frame=i_map2.coordinate_frame),
                      top_right=SkyCoord(*top_right, frame=i_map2.coordinate_frame))

fig = plt.figure()
ax = fig.add_subplot(projection=i_map)
i_map.plot(axes=ax, clip_interval=(1, 99.9) * u.percent, zorder=0, autoalign=True, cmap='RdYlBu')
i_map1.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=1, autoalign=True, cmap='RdYlBu')
i_map2.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=2, autoalign=True, cmap='RdYlBu')
#i_map3.plot(axes=ax, clip_interval=(1, 99.97) * u.percent, alpha=0.7, zorder=3, autoalign=True, cmap='RdYlBu')

aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Width Composite Map")
plt.show()
