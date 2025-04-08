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
#11
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_012613/EIS_fit_fe_12_195.asdf") as af10:
    intensity_map10 = af10["int_map"].data
    i_map10 = Map(intensity_map10, af10["int_map"].meta)
#12
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_163612/EIS_fit_fe_12_195.asdf") as af11:
    intensity_map11 = af11["int_map"].data
    i_map11 = Map(intensity_map11, af11["int_map"].meta)
#13
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_030043/EIS_fit_fe_12_195.asdf") as af12:
    intensity_map12 = af12["int_map"].data
    i_map12 = Map(intensity_map12, af12["int_map"].meta)
#14
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_043942/EIS_fit_fe_12_195.asdf") as af13:
    intensity_map13 = af13["int_map"].data
    i_map13 = Map(intensity_map13, af13["int_map"].meta)
#15
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_061712/EIS_fit_fe_12_195.asdf") as af14:
    intensity_map14 = af14["int_map"].data
    i_map14 = Map(intensity_map10, af14["int_map"].meta)
#16
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_075442/EIS_fit_fe_12_195.asdf") as af15:
    intensity_map15 = af15["int_map"].data
    i_map15 = Map(intensity_map15, af15["int_map"].meta)
#17
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_093012/EIS_fit_fe_12_195.asdf") as af10:
    intensity_map16 = af16["int_map"].data
    i_map16 = Map(intensity_map16, af16["int_map"].meta)
#18
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_105442/EIS_fit_fe_12_195.asdf") as af17:
    intensity_map17 = af17["int_map"].data
    i_map17 = Map(intensity_map17, af17["int_map"].meta)
#19
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_120542/EIS_fit_fe_12_195.asdf") as af18:
    intensity_map18 = af18["int_map"].data
    i_map18 = Map(intensity_map18, af18["int_map"].meta)
#20
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_131642/EIS_fit_fe_12_195.asdf") as af19:
    intensity_map19 = af19["int_map"].data
    i_map19 = Map(intensity_map19, af19["int_map"].meta)
#21
with asdf.open("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_144842/EIS_fit_fe_12_195.asdf") as af20:
    intensity_map20 = af20["int_map"].data
    i_map20 = Map(intensity_map20, af20["int_map"].meta)

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
i_smap11 = i_map11.submap(SkyCoord(*bottom_left, frame=i_map11.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map11.coordinate_frame))
i_smap12 = i_map12.submap(SkyCoord(*bottom_left, frame=i_map12.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map12.coordinate_frame))
i_smap13 = i_map13.submap(SkyCoord(*bottom_left, frame=i_map13.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map13.coordinate_frame))
i_smap14 = i_map14.submap(SkyCoord(*bottom_left, frame=i_map14.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map14.coordinate_frame))
i_smap15 = i_map15.submap(SkyCoord(*bottom_left, frame=i_map15.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map15.coordinate_frame))
i_smap16 = i_map16.submap(SkyCoord(*bottom_left, frame=i_map16.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map16.coordinate_frame))
i_smap17 = i_map17.submap(SkyCoord(*bottom_left, frame=i_map17.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map17.coordinate_frame))
i_smap18 = i_map18.submap(SkyCoord(*bottom_left, frame=i_map18.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map18.coordinate_frame))
i_smap19 = i_map19.submap(SkyCoord(*bottom_left, frame=i_map19.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map19.coordinate_frame))
i_smap20 = i_map20.submap(SkyCoord(*bottom_left, frame=i_map20.coordinate_frame),
                        top_right=SkyCoord(*top_right, frame=i_map20.coordinate_frame))

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

aspect_ratio = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
ax.set_aspect(aspect_ratio)

ax.set_title("EIS Iron XII Composite Map")
plt.show()



