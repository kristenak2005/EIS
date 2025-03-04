import matplotlib.pyplot as plt
import matplot.image as mpimg

img_path = '/mnt/scratch/data/spruksk2/python_output/EIS_work/20151018_102719/plots/EIS_composition_20151018_102719_CaAr.png'
img = mpimg.imread(img_path)
imgplot plt.imshow(img)
plt.show()
