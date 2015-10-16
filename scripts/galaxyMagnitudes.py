import numpy as np
from scipy.optimize import curve_fit
import dataIO
import matplotlib.pyplot as plt
from matplotlib import rc


# change typesetting for plots
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


#define two gaussians
def twoGauss(x, a1, c1, s1, a2, c2, s2):
    if s1 ==0 or s2 == 0:
        return 99999.9
    else:
        return a1*np.exp(-(x-c1)**2/(2*s1**2))+a2*np.exp(-(x-c2)**2/(2*s2**2))



# create a histogram of all of the galaxy magnitudes
def mag_histogram():
    # get the dictionary of galaxies
    galaxies = dataIO.galaxy_init()

    mags = np.array([])
    brightzs = np.array([])
    dimzs = np.array([])
    # loop through all the galaxies and add the magnitudes to a list
    for key in galaxies:
        if galaxies[key].R > 0:
            mags = np.append(mags, galaxies[key].R)
            if galaxies[key].R >= 27:
                dimzs = np.append(dimzs, galaxies[key].z)
            else:
                brightzs = np.append(brightzs, galaxies[key].z)


    # fit both dim and bright galaxies to a double gauss
    brightHist, bin_edges = np.histogram(brightzs, range=(3.03, 3.15), bins=31)
    dimHist, bin_edges = np.histogram(dimzs, range=(3.03, 3.15), bins=31)
    # find the centers of the bins
    bin_centers = [0.5*(bin_edges[ii]+bin_edges[ii+1]) for ii in range(len(bin_edges)-1)]

    fitParamBright, fitCovBright = curve_fit(twoGauss, bin_centers, brightHist)
    fitParamDim, fitCovDim = curve_fit(twoGauss, bin_centers, dimHist)

    # create a new x array
    x = np.linspace(3.03, 3.15, 500)

#    plt.plot(x, twoGauss(x, *fitParamBright), 'r', linewidth=2)
#    plt.plot(x, twoGauss(x, *fitParamDim), 'b', linewidth=2)

    print()
    # find which peak is the red and which peak is the blue peak
    # for the bright objects
    if fitParamBright[1] > fitParamBright[4]:
        blue = 1
        red = 0
    else:
        blue = 0
        red = 1

    # first look at the blue peak
    print("Bright Blue Peak Parameters:")
    print("Center: z={:.4f}+/-{:.4f}".format(fitParamBright[1+3*blue],np.sqrt(np.abs(fitCovBright[1+3*blue][1+3*blue])) ))
    
    print("Amplitude: z={:.4f}+/-{:.4f}".format(fitParamBright[0+3*blue],np.sqrt(np.abs(fitCovBright[0+3*blue][0+3*blue])) ))
    print("Width: z={:.4f}+/-{:.4f}".format(fitParamBright[2+3*blue],np.sqrt(np.abs(fitCovBright[2+3*blue][2+3*blue])) ))

    print()
    print("Bright Red Peak Parameters:")
    print("Center: z={:.4f}+/-{:.4f}".format(fitParamBright[1+3*red],np.sqrt(np.abs(fitCovBright[1+3*red][1+3*red])) ))
    
    print("Amplitude: z={:.4f}+/-{:.4f}".format(fitParamBright[0+3*red],np.sqrt(np.abs(fitCovBright[0+3*red][0+3*red])) ))
    print("Width: z={:.4f}+/-{:.4f}".format(fitParamBright[2+3*red],np.sqrt(np.abs(fitCovBright[2+3*red][2+3*red])) ))


    print()
    # find which peak is the red and which peak is the blue peak
    # for the bright objects
    if fitParamDim[1] > fitParamDim[4]:
        blue = 1
        red = 0
    else:
        blue = 0
        red = 1

    # first look at the blue peak
    print("Dim Blue Peak Parameters:")
    print("Center: z={:.4f}+/-{:.4f}".format(fitParamDim[1+3*blue],np.sqrt(np.abs(fitCovDim[1+3*blue][1+3*blue])) ))
    
    print("Amplitude: z={:.4f}+/-{:.4f}".format(fitParamDim[0+3*blue],np.sqrt(np.abs(fitCovDim[0+3*blue][0+3*blue])) ))
    print("Width: z={:.4f}+/-{:.4f}".format(fitParamDim[2+3*blue],np.sqrt(np.abs(fitCovDim[2+3*blue][2+3*blue])) ))

    print()
    print("Dim Red Peak Parameters:")
    print("Center: z={:.4f}+/-{:.4f}".format(fitParamDim[1+3*red],np.sqrt(np.abs(fitCovDim[1+3*red][1+3*red])) ))
    
    print("Amplitude: z={:.4f}+/-{:.4f}".format(fitParamDim[0+3*red],np.sqrt(np.abs(fitCovDim[0+3*red][0+3*red])) ))
    print("Width: z={:.4f}+/-{:.4f}".format(fitParamDim[2+3*red],np.sqrt(np.abs(fitCovDim[2+3*red][2+3*red])) ))








    print()


    plt.hist(brightzs, range=(3.03, 3.15), histtype="stepfilled", linewidth=2,  bins=30, label="$M_R \ < \ 27$", color='red')
    plt.hist(dimzs, range=(3.03, 3.15), histtype="stepfilled", linewidth=2,  bins=30, label="$M_R \ \geq \ 27$", alpha=0.6)
    plt.legend(loc="upper left")
    plt.xlabel("$z$", fontsize=16)
    plt.ylabel("$N$", fontsize=16)
    plt.xlim([3.03, 3.12])
    plt.savefig("histogram_Rseparate.png", dpi=400)
#    plt.show()








if __name__=="__main__":
    mag_histogram()
