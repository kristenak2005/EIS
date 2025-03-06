import matplotlib.pyplot as plt
import asdf

with asdf.open("EIS_fit_fe_12_195.asdf") as af:
    intensity_map = af["int_map"].data
    plt.imshow(intensity_map. cmap="RdYlBu", origin="lower")
    plt.title("EIS fit FE XII 195")
    plt.show()


