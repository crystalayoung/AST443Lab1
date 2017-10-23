# This code is used for plotting data used to analyze flux

import numpy as np
import matplotlib.pyplot as plt

#sys.path.append('Users/crystalyoung/Desktop/AST_443_Lab_0')

def plot_flux(text):
    text = np.loadtxt(text)
    x = text[:,0]
    y = text[:,1]
    plt.plot(x,y)
    plt.xlabel('Wavelength (Angstroms)')
    plt.title('Spectrum of Arc Lamp')
    plt.show()