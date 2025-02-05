import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u

# Load your EIS FITS file
eis_map = sunpy.map.Map("eis_20151018_102719_sis.fits")

# Plot the image with clipping (optional)
eis_map.peek(clip_interval=(1, 99.99) * u.percent)

# Create a dictionary (ASDF tree) and save the file
tree = {'eis_map': eis_map}
with asdf.AsdfFile(tree) as asdf_file:
    asdf_file.write_to("eis_map.asdf")

print("ASDF file saved successfully!")
