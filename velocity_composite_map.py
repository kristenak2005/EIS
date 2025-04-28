import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import os
from astropy.visualization import ImageNormalize, LinearStretch

eis_evts = ['20151018_102719','20151018_113839','20151018_124939','20151018_140939','20151018_160113','20151018_173743',
            '20151018_191443','20151018_205143','20151018_222743','20151018_235413','20151019_012613','20151019_030043',
            '20151019_043942','20151019_061712','20151019_075442','20151019_093012','20151019_105442','20151019_120542',
            '20151019_131642','20151019_144842','20151019_163612','20151019_181312','20151019_194942','20151019_212642',
            '20151019_230142']

# Set min/max velocity (adjust as needed)
min_value = -50  # km/s
max_value = 50   # km/s

# Open first velocity map to define figure setup
v_map = Map(f'/mnt/scratch/data/spruksk2/python_output/EIS_work/{eis_evts[0]}/save_files/velocity_files/eis_{eis_evts[0]}_velocity_fe_12_195.fits')

fig = plt.figure()
ax = fig.add_subplot(projection=v_map)

# Define FOV in world coordinates
xlims = [-1000, 1000] * u.arcsec
ylims = [-1000, 1000] * u.arcsec
world_coords = SkyCoord(Tx=xlims, Ty=ylims, frame=v_map.coordinate_frame)
pixel_coords_x, pixel_coords_y = v_map.wcs.world_to_pixel(world_coords)

# Set normalization
norm = ImageNormalize(vmin=min_value, vmax=max_value, stretch=LinearStretch())

# Start layering velocity maps
z_val = 0
for evt in eis_evts:
    v_map = Map(f'/mnt/scratch/data/spruksk2/python_output/EIS_work/{evt}/save_files/velocity_files/eis_{evt}_velocity_fe_12_195.fits')
    aspect = v_map.meta['cdelt2'] / v_map.meta['cdelt1']
    v_map.plot(axes=ax, norm=norm, alpha=0.7, zorder=z_val, autoalign=True, cmap='RdBu_r', aspect=aspect)
    z_val += 1

ax.set_xlim(pixel_coords_x)
ax.set_ylim(pixel_coords_y)

# Save output
output_filename = "/mnt/scratch/data/spruksk2/python_output/composite_velocity_map.png"
plt.savefig(output_filename, bbox_inches='tight')

ax.set_title("EIS Velocity Composite Map")
plt.show()
