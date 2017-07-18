from hdf_adam import *
import matplotlib.pyplot as plt

filename = raw_input("Filename: " )

h5_data = read_hdf(filename)

plt.figure()
plt.imshow(h5_data.data,aspect='auto',extent=[h5_data.axes[0].axis_min,
	h5_data.axes[0].axis_max,h5_data.axes[1].axis_min,h5_data.axes[1].axis_max],norm=matplotlib.colors.SymLogNorm(10**-6),cmap='Spectral',origin='lower')
plt.colorbar()
plt.xlabel("%s \n" % (h5_data.axes[0].attributes['LONG_NAME'][0]))
plt.ylabel("%s \n" % (h5_data.axes[1].attributes['LONG_NAME'][0]))
plt.tight_layout()
plt.show()

 
