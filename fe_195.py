import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import asdf

#1
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/EIS_fit_fe_12_195.asdf") as af:
    intensity_map = af["int_map"].data
    i_map = Map(intensity_map, af["int_map"].meta)
#2
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/EIS_fit_fe_12_195.asdf") as af1:
    intensity_map1 = af1["int_map"].data
    i_map1 = Map(intensity_map1, af1["int_map"].meta)
#3
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/EIS_fit_fe_12_195.asdf") as af2:
    intensity_map2 = af2["int_map"].data
    i_map2 = Map(intensity_map2, af2["int_map"].meta)
#4
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_140939/EIS_fit_fe_12_195.asdf") as af3:
    intensity_map3 = af3["int_map"].data
    i_map3 = Map(intensity_map3, af3["int_map"].meta)
#5
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_160113/EIS_fit_fe_12_195.asdf") as af4:
    intensity_map4 = af4["int_map"].data
    i_map4 = Map(intensity_map4, af4["int_map"].meta)
#6
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_173743/EIS_fit_fe_12_195.asdf") as af5:
    intensity_map5 = af5["int_map"].data
    i_map5 = Map(intensity_map5, af5["int_map"].meta)
#7
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/EIS_fit_fe_12_195.asdf") as af6:
    intensity_map6 = af6["int_map"].data
    i_map6 = Map(intensity_map6, af6["int_map"].meta)
#8
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/EIS_fit_fe_12_195.asdf") as af7:
    intensity_map7 = af7["int_map"].data
    i_map7 = Map(intensity_map7, af7["int_map"].meta)
#9
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_222743/EIS_fit_fe_12_195.asdf") as af8:
    intensity_map8 = af8["int_map"].data
    i_map8 = Map(intensity_map8, af8["int_map"].meta)
#10
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_235413/EIS_fit_fe_12_195.asdf") as af9:
    intensity_map9 = af9["int_map"].data
    i_map9 = Map(intensity_map9, af9["int_map"].meta)



# Define the region of interest using SkyCoord
bottom_left = [200, -800] * u.arcsec
top_right = [1000, -200] * u.arcsec

# Create submaps for both maps
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
i_smap7 = i_map7.submap(SkyCoord(*bottom_left, frame=i_map7.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map7.coordinate_frame))
i_smap8 = i_map8.submap(SkyCoord(*bottom_left, frame=i_map8.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map8.coordinate_frame))
i_smap9 = i_map9.submap(SkyCoord(*bottom_left, frame=i_map9.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map9.coordinate_frame))
i_smap10 = i_map10.submap(SkyCoord(*bottom_left, frame=i_map10.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map10.coordinate_frame))

# Create a figure and axis for plotting
fig = plt.figure()
ax = fig.add_subplot(projection=i_map)

i_map.plot(axes=ax, clip_interval=(1, 99.9) * u.percent, zorder=0, autoalign=True, cmap='RdYlBu')
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


aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Iron XII Composite Map")
plt.show()



