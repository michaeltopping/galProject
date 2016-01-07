import numpy as np
import matplotlib.pyplot as plt
import pyfits
from matplotlib import rc




# change the fonts for the plotting windows
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})
rc('text', usetex=True)






# a function to plot a single spectrum
# input:
#  obj - the name of the object that is plotted
#  path - the path to the spectrum data
#  filename - the filename of the spectrum data file
def single_plot(obj, path, filename):
    print("Plotting object: {}".format(obj))


    # read in the spectrum data
    spec, header = pyfits.getdata(path+filename, 0, header=True)



    # compute the wavelength solution
    if 'CRVAL1' in header:
        lambdaMin = float(header['CRVAL1'])
        dlambda = float(header['CD1_1'])
        Nlambda = float(header['NAXIS1'])

    # if the CRVAL1 parameters are not in the header file, we have to use a different
    #  wavelength solution
    else:
        # here we will loop through parameters in the WAT2_002 header key and look 
        #  and look for the first one that is large enough to be a wavelength limit
        #  then we can find the limit, the dlambda from the other parameters
        foundMin = False
        index = 0
        while not(foundMin):
            if 'WAT2_002' in header:
                print("One element: ", float(header['WAT2_002'].split()[index]))
                if float(header['WAT2_002'].split()[index]) > 2000.0:
                    founrMin = True
            
                

            else:
                pass


            lambdaMin = float(header['WAT2_002'].split()[index])
            dlambda = float(header['WAT2_002'].split()[index+1])
            Nlambda = np.shape(spec)[0]




    # create the wavelengths array
    wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)



    # do the plotting stuff here
    plot_hist(wavelengths, spec)
    plt.xlabel('$\lambda$')
    plt.xlim([4800, 5100])
    
    plt.show()




# plot multiple spectra on one figure
def multi_plot(objs, paths, filenames, axs):
    for obj, path, filename, ax in zip(objs, paths, filenames, axs):
        print("Plotting object: {}".format(obj))


        # read in the spectrum data
        spec, header = pyfits.getdata(path+filename, 0, header=True)



        # compute the wavelength solution
        if 'CRVAL1' in header:
            lambdaMin = float(header['CRVAL1'])
            dlambda = float(header['CD1_1'])
            Nlambda = float(header['NAXIS1'])

        # if the CRVAL1 parameters are not in the header file, we have to use a different
        #  wavelength solution
        else:
            # here we will loop through parameters in the WAT2_002 header key and look 
            #  and look for the first one that is large enough to be a wavelength limit
            #  then we can find the limit, the dlambda from the other parameters
            foundMin = False
            index = 0
            while not(foundMin):
                if 'WAT2_002' in header:
                    print("One element: ", float(header['WAT2_002'].split()[index]))
                    if float(header['WAT2_002'].split()[index]) > 2000.0:
                        foundMin = True
                
                    

                else:
                    pass


                lambdaMin = float(header['WAT2_002'].split()[index])
                dlambda = float(header['WAT2_002'].split()[index+1])
                Nlambda = np.shape(spec)[0]




        # create the wavelengths array
        wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)




        #set limits on the wavelength that we are looking at
        minWavelength = 4920
        maxWavelength = 5020
        minWavelengthIndex = np.where(abs(wavelengths-minWavelength) ==
                               min(abs(wavelengths-minWavelength)))[0]
        maxWavelengthIndex = np.where(abs(wavelengths-maxWavelength) ==
                               min(abs(wavelengths-maxWavelength)))[0]

        print(minWavelengthIndex, maxWavelengthIndex)
        wavelengths = wavelengths[minWavelengthIndex:maxWavelengthIndex]
        spec = spec[minWavelengthIndex:maxWavelengthIndex]



        # do the plotting stuff here
        plot_hist(wavelengths, spec, ax, obj)
        plt.xlabel(r'Wavelength $[\rm \AA]$')
        plt.xlim([4920, 5020])
        
        # put text information here, like object name
        ax.text(0.03, 0.85, obj,  transform=ax.transAxes, fontsize=16) 






# plot a spectrum as a histogram
def plot_hist(wavelengths, spec, axes, obj):
    
    # loop through each wavelength
    for ii in range(1,wavelengths.size-1):
        axes.plot([wavelengths[ii-1], wavelengths[ii-1]], [spec[ii-1], spec[ii]], 'k', linewidth=2)
        axes.plot([wavelengths[ii-1], wavelengths[ii]], [spec[ii], spec[ii]], 'k', linewidth=2)

    # get rid of the y axis ticks
    yticks = axes.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)
    

    axes.set_yticklabels([])

    # plot the ylabel
    axes.set_ylabel('Flux')
    

# main function where the plotting will be called
if __name__ =="__main__":
#    single_plot("C11","../shapley2006_spec/","ssa22.C11.br.fits")

    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)


    multi_plot(["C11", r"laec4\text{_}0721", "NB2547"], ["../shapley2006_spec/", "../spec/", "../spec/"], ["ssa22.C11.br.fits", "ssaly_8.b.laec4_0721.msdfc_v.fits", "ssa22_nb1.b.NB2547.msdfc_v.fits"], [ax1, ax2, ax3])


    # make the plots have no space inbetween them
    f.subplots_adjust(hspace=0)


    # save the figure
    plt.savefig("multispec.png", dpi=400)

    # show the plot 
#    plt.show()
