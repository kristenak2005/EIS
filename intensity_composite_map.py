import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import os
import asdf
from astropy.visualization import ImageNormalize, LinearStretch
eis_evts = ['20151018_102719','20151018_113839','20151018_124939','20151018_140939','20151018_160113','20151018_173743',
            '20151018_191443','20151018_205143','20151018_222743','20151018_235413','20151019_012613','20151019_030043',
            '20151019_043942','20151019_061712','20151019_075442','20151019_093012','20151019_105442','20151019_120542',
            '20151019_131642','20151019_144842','20151019_163612','20151019_181312','20151019_194942','20151019_212642',
            '20151019_230142']
min_value = 0
max_value = 2000

arrs = asdf.open('/mnt/scratch/data/spruksk2/python_output/EIS_work/'+eis_evts[0]+'/EIS_fit_fe_12_195.asdf')
i_map = arrs['int_map']

fig = plt.figure()
ax = fig.add_subplot(projection=i_map)

xlims = [-1000,1000]*u.arcsec
ylims = [-1000,1000]*u.arcsec
world_coords = SkyCoord(Tx=xlims, Ty=ylims, frame=i_map.coordinate_frame)
pixel_coords_x, pixel_coords_y = i_map.wcs.world_to_pixel(world_coords)

#norm = colors.Normalize(vmin = min_value, vmax = max_value)
norm = ImageNormalize(vmin=min_value, vmax=max_value, stretch=LinearStretch())

z_val = 0

for evt in eis_evts:
    arrs = asdf.open('/mnt/scratch/data/spruksk2/python_output/EIS_work/'+evt+'/EIS_fit_fe_12_195.asdf')
    i_map = arrs['int_map']
    aspect = i_map.meta['cdelt2'] / i_map.meta['cdelt1']
    i_map.plot(axes=ax, norm=norm, alpha=0.7, zorder=z_val, autoalign=True, cmap='OrRd', aspect=aspect)
    z_val += 1


bottom_left_active = SkyCoord((300) * u.arcsec, (-550) * u.arcsec, frame=i_map.coordinate_frame)
top_right_active = SkyCoord((700) * u.arcsec, (-250) * u.arcsec, frame=i_map.coordinate_frame)

# Quiet region centered at (-450, 500)
bottom_left_quiet = SkyCoord((100) * u.arcsec, (0) * u.arcsec, frame=i_map.coordinate_frame)
top_right_quiet = SkyCoord((500) * u.arcsec, (300) * u.arcsec, frame=i_map.coordinate_frame)

# Draw quadrangles for active and quiet regions
i_map.draw_quadrangle(bottom_left_active, axes=ax, top_right=top_right_active, edgecolor="red")
i_map.draw_quadrangle(bottom_left_quiet, axes=ax, top_right=top_right_quiet, edgecolor="black")

ax.set_xlim(pixel_coords_x)
ax.set_ylim(pixel_coords_y)

ax.set_xlabel('Solar X [arcsec]', fontsize=12)
ax.set_ylabel('Solar Y [arcsec]', fontsize=12)
ax.set_title("EIS Intensity (Fe XII 195 Å) Composite Map", fontsize=14)
ax.tick_params(labelsize=10)


#save_dir = "mnt/scratch/data/spruksk2/python_output/EIS_work/"
#save_dir = "/Users/dml/python_output/"
#output_filename = os.path.join(save_dir, "composite_intensity_map.png")
#plt.savefig("/mnt/scratch/data/spruksk2/python_output/composite_intensity_map.png", bbox_inches='tight')
                              
ax.set_title("EIS Intensity Composite Map")
plt.show()
