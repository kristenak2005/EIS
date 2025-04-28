import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import os
import asdf
from astropy.visualization import ImageNormalize, LinearStretch

#i_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/intensity_files/eis_20151018_102719_intensity_fe_12_195.fits")   
#i_map1 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_113839/save_files/intensity_files/eis_20151018_113839_intensity_fe_12_195.fits") 
#i_map2 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_124939/save_files/intensity_files/eis_20151018_124939_intensity_fe_12_195.fits")
#i_map3 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_140939/save_files/intensity_files/eis_20151018_140939_intensity_fe_12_195.fits")
#i_map4 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_160113/save_files/intensity_files/eis_20151018_160113_intensity_fe_12_195.fits")
#i_map5 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_173743/save_files/intensity_files/eis_20151018_173743_intensity_fe_12_195.fits")
#i_map6 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_191443/save_files/intensity_files/eis_20151018_191443_intensity_fe_12_195.fits")
#i_map7 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_205143/save_files/intensity_files/eis_20151018_205143_intensity_fe_12_195.fits")
#i_map8 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_222743/save_files/intensity_files/eis_20151018_222743_intensity_fe_12_195.fits")
#i_map9 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_235413/save_files/intensity_files/eis_20151018_235413_intensity_fe_12_195.fits")
#i_map10 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_012613/save_files/intensity_files/eis_20151019_012613_intensity_fe_12_195.fits")
#i_map11 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_030043/save_files/intensity_files/eis_20151019_030043_intensity_fe_12_195.fits")
#i_map12 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_043942/save_files/intensity_files/eis_20151019_043942_intensity_fe_12_195.fits")
#i_map13 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_061712/save_files/intensity_files/eis_20151019_061712_intensity_fe_12_195.fits")
#i_map14 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_075442/save_files/intensity_files/eis_20151019_075442_intensity_fe_12_195.fits")
#i_map15 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_093012/save_files/intensity_files/eis_20151019_093012_intensity_fe_12_195.fits")
#i_map16 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_105442/save_files/intensity_files/eis_20151019_105442_intensity_fe_12_195.fits")
#i_map17 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_120542/save_files/intensity_files/eis_20151019_120542_intensity_fe_12_195.fits")
#i_map18 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_131642/save_files/intensity_files/eis_20151019_131642_intensity_fe_12_195.fits")
#i_map19 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_144842/save_files/intensity_files/eis_20151019_144842_intensity_fe_12_195.fits")
#i_map20 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_163612/save_files/intensity_files/eis_20151019_163612_intensity_fe_12_195.fits")
#i_map21 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_181312/save_files/intensity_files/eis_20151019_181312_intensity_fe_12_195.fits")
#i_map22 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_194942/save_files/intensity_files/eis_20151019_194942_intensity_fe_12_195.fits")
#i_map23 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_212642/save_files/intensity_files/eis_20151019_212642_intensity_fe_12_195.fits")
#i_map24 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151019_230142/save_files/intensity_files/eis_20151019_230142_intensity_fe_12_195.fits")
#i_map25 = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151020_002612/save_files/intensity_files/eis_20151020_002612_intensity_fe_12_195.fits")
eis_evts = ['20151018_102719','20151018_113839','20151018_124939','20151018_140939','20151018_160113','20151018_173743',
            '20151018_191443','20151018_205143','20151018_222743','20151018_235413','20151019_012613','20151019_030043',
            '20151019_043942','20151019_061712','20151019_075442','20151019_093012','20151019_105442','20151019_120542',
            '20151019_131642','20151019_144842','20151019_163612','20151019_181312','20151019_194942','20151019_212642',
            '20151019_230142']


#v_map = Map("/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_velocity_fe_12_195.fits")
#w_map = Map('/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/save_files/velocity_files/eis_20151018_102719_width_fe_12_195.fits')


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

ax.set_xlim(pixel_coords_x)
ax.set_ylim(pixel_coords_y)


#save_dir = "mnt/scratch/data/spruksk2/python_output/EIS_work/"
#save_dir = "/Users/dml/python_output/"
#output_filename = os.path.join(save_dir, "composite_intensity_map.png")
plt.savefig("/mnt/scratch/data/spruksk2/python_output/composite_intensity_map.png", bbox_inches='tight')
                              
ax.set_title("EIS Intensity Composite Map")
plt.show()
