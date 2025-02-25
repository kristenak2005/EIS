import matplotlib.pyplot as plt
import asdf

with asdf.open("EIS_fit_fe_12_195.asdf") as af:
    intensity_map = af["int_map"].data
    plt.imshow(intensity_map. cmap="inferno", origin="lower")
    plt.title("EIS fit FE XII 195")
    plt.show()




# Plot with better scaling
plt.figure(figsize=(8, 6))  # Adjust figure size
plt.imshow(intensity_map, cmap="inferno", origin="lower", aspect="equal")  # ✅ Maintain aspect ratio
plt.colorbar(label="Intensity")  # Add a colorbar
plt.title("EIS Fit Fe XII 195")
plt.xlabel("X Pixel")
plt.ylabel("Y Pixel")
plt.show()
