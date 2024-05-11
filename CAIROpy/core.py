# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


class CAIRO(object):
    def __init__(self, cube, hdr=None):
        super(CAIRO, self).__init__()
        self.cube = cube
        self.hdr = hdr if hdr is not None else None

        
    def rms_map(self, rms_map, filename=None):
        if not filename :
            print("Generate rms_map.dat file")
        else: print("Generate " + filename + " file")

        filename = filename or "myrms_map.dat"
        
        with open(filename,'w+') as f:
            f.write('{:d}\t{:d}\n'.format(self.cube.shape[1], self.cube.shape[2]))
            for j in range(self.cube.shape[2]):
                for k in range(self.cube.shape[1]):
                    line = '{:d}\t{:d}\t{:0.16f}\n'.format(k, j, rms_map[k,j])
                    f.write(line)
        

    def gen_parameters(self, filename_parameters=None, filename=None, fileout="result.dat", timeout="timestep.dat", filename_noise="",
                       filename_init_spec="", filename_lsf="", n_gauss=6, lambda_amp=1, lambda_mu=1, lambda_sig=1,
                       lambda_sig_corr_narrow=1, lambda_sig_corr_broad=0, lambda_mu_corr_narrow=1, lambda_mu_corr_broad=0,
                       lambda_var_sig=1, lambda_r=1, lb_sig=1, ib_sig=20, ub_sig=60, lb_amp=0.01, ub_amp=0,
                       maxiter=800, iprint = -1, descent=".false."):
        if not filename :
            print("Need an input filename")
            sys.exit()
            
        if not filename_parameters :
            print("Generate parameters.txt file")
        else: print("Generate " + filename_parameters + " file")

        filename_parameters = filename_parameters or "parameters.txt"
        
        input_file = open(filename_parameters, 'w')
        input_file.write("&user_parameters"+'\n')
        input_file.write("    filename =  "+repr(filename)+'\n')
        input_file.write("    ,fileout =  "+repr(fileout)+'\n')
        input_file.write("    ,timeout =  "+repr(timeout)+'\n')
        input_file.write("    ,filename_noise =  "+repr(filename_noise)+'\n')
        input_file.write("    ,filename_init_spec =  "+repr(filename_init_spec)+'\n')
        input_file.write("    ,filename_lsf =  "+repr(filename_lsf)+'\n')
        input_file.write("    ,n_gauss =  "+repr(n_gauss)+'\n')
        input_file.write("    ,lambda_amp =  "+repr(lambda_amp)+'d0'+'\n')
        input_file.write("    ,lambda_mu =  "+repr(lambda_mu)+'d0'+'\n')
        input_file.write("    ,lambda_sig =  "+repr(lambda_sig)+'d0'+'\n')
        input_file.write("    ,lambda_sig_corr_narrow =  "+repr(lambda_sig_corr_narrow)+'d0'+'\n')
        input_file.write("    ,lambda_sig_corr_broad =  "+repr(lambda_sig_corr_broad)+'d0'+'\n')
        input_file.write("    ,lambda_mu_corr_narrow =  "+repr(lambda_mu_corr_narrow)+'d0'+'\n')
        input_file.write("    ,lambda_mu_corr_broad =  "+repr(lambda_mu_corr_broad)+'d0'+'\n')
        input_file.write("    ,lambda_var_sig =  "+repr(lambda_var_sig)+'d0'+'\n')
        input_file.write("    ,lambda_r =  "+repr(lambda_r)+'d0'+'\n')
        input_file.write("    ,lb_sig =  "+repr(lb_sig)+'d0'+'\n')
        input_file.write("    ,ib_sig =  "+repr(ib_sig)+'d0'+'\n')
        input_file.write("    ,ub_sig =  "+repr(ub_sig)+'d0'+'\n')
        input_file.write("    ,lb_amp =  "+repr(lb_amp)+'d0'+'\n')
        input_file.write("    ,ub_amp =  "+repr(ub_amp)+'d0'+'\n')
        input_file.write("    ,maxiter =  "+repr(maxiter)+'\n')
        input_file.write("    ,iprint =  "+repr(iprint)+'\n')
        input_file.write("    ,descent =  "+descent+'\n')
        input_file.write("    /"+'\n')
        input_file.close()


    def run(self, filename=None, nohup=False):
        if not filename: 
            print("Need the input filename parameters to run CAIRO")
            sys.exit()
        if nohup == False:
            os.system("CAIRO " + filename)
        else:
            os.system("nohup CAIRO " + filename + "&")


    def physical_gaussian(self, gaussian):
        n_gauss = gaussian.shape[0]/3
        output = np.zeros(gaussian.shape)
        if self.hdr is not None :
            output[0::3] = gaussian[0::3]
            output[1::3] = self.mean2vel(self.hdr["CRVAL3"]*1.e-3, self.hdr["CDELT3"]*1.e-3, self.hdr["CRPIX3"]-1, gaussian[1::3])
            output[2::3] = gaussian[2::3] * np.abs(self.hdr["CDELT3"])*1.e-3
            return output
        else:
            print("Missing header")
            return 0        


    def gauss(self, x, a, mu, sig):
        return a * np.exp(-((x - mu)**2)/(2. * sig**2))


    def gauss_2D(self, xs, a, mu, sig):
        return [a * np.exp(-((x - mu)**2)/(2. * sig**2)) for x in xs]

    
    def mean2vel(self, CRVAL, CDELT, CRPIX, mean):
        return [(CRVAL + CDELT * (mean[i] - CRPIX)) for i in range(len(mean))]            


    def plot_spect(self, gaussian, idy=0, idx=0):
        n_gauss = gaussian.shape[0]/3
        x = np.arange(self.cube.shape[0])

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        if self.hdr is not None :
            if not self.hdr["CRVAL3"] : print("Missing CRVAL3 keyword")
            v = self.mean2vel(self.hdr["CRVAL3"]*1.e-3, self.hdr["CDELT3"]*1.e-3, self.hdr["CRPIX3"]-1, x)
            if v[0] > v[1] : v = v[::-1]
            ax.step(v, self.cube[:,idy,idx], color='cornflowerblue')
            tot = np.zeros(self.cube.shape[0])
            for i in np.arange(n_gauss):
                spectrum = self.gauss(x, gaussian[int(0+(3*i)),idy,idx], gaussian[int(1+(3*i)),idy,idx], gaussian[int(2+(3*i)),idy,idx])
                tot += spectrum
                ax.plot(v, spectrum, color="k")
            ax.plot(v, tot, color="r") 
            ax.set_ylabel(r'T [k]')
            ax.set_xlabel(r'v [km s$^{-1}$]')
        else:
            ax.step(x, self.cube[:,idy,idx], color='cornflowerblue')
            tot = np.zeros(self.cube.shape[0])
            for i in np.arange(n_gauss):
                spectrum = self.gauss(x, gaussian[int(0+(3*i)),idy,idx], gaussian[int(1+(3*i)),idy,idx], gaussian[int(2+(3*i)),idy,idx])
                tot += spectrum
                ax.plot(x, spectrum, color="k")
            ax.plot(x, tot, color="r") 
            ax.set_ylabel(r'T [k]')
            ax.set_xlabel(r'idx [pixel unit]')
             
        return 0 


    def return_result_cube(self, gaussian=None, ampfield=None, pixfield=None, sigfield=None):
        if gaussian is not None:
            result = np.zeros(self.cube.shape)
            n_gauss = gaussian.shape[0]/3
            for i in np.arange(n_gauss) :
                result += self.gauss_2D(np.arange(self.cube.shape[0]), gaussian[int(0+(3*i))], gaussian[int(1+(3*i))], gaussian[int(2+(3*i))])
            return result
        elif ampfield is not None and pixfield is not None and sigfield is not None:
            result = np.zeros(self.cube.shape)
            n_gauss = ampfield.shape[0]
            for i in np.arange(n_gauss) :
                result += self.gauss_2D(np.arange(self.cube.shape[0]), ampfield[i], pixfield[i], sigfield[i])
            return result
        else: print("error : 1 or 3 arguments needed")
                
        
if __name__ == '__main__':    
    #Load data
    path="/home/amarchal/CAIRO/data/"
    filename = "N5236_cropped_subcube.fits"
    hdu = fits.open(path+filename)
    hdr = hdu[0].header
    cube = hdu[0].data

    #Call CAIROpy
    core = CAIRO(cube, hdr=hdr)

    filename_parameters = path + "N5236_cropped_subcube_parameters_run_0.txt"
    filename = path + "N5236_cropped_subcube.fits"
    fileout = "!" + path + "N5236_cropped_subcube_gauss_run_0.fits"
    filename_noise = path + "N5236_cropped_subcube_rms.fits"
    filename_init_spec = path + "N5236_input.csv"
    filename_lsf = path + "N5236_lsf.csv"
    n_gauss = 6
    lambda_amp = 10.
    lambda_mu = 10.
    lambda_sig = 10.
    lambda_sig_corr_narrow = 10.
    lambda_sig_corr_broad = 0.
    lambda_mu_corr_narrow = 10.
    lambda_mu_corr_broad = 0.
    lambda_var_sig = 1.
    lambda_r = 1.
    lb_sig = 1.
    ib_sig = 4.
    ub_sig = 10.
    lb_amp = 0.01
    ub_amp = 0.
    maxiter = 800
    iprint = 25
    descent = ".false."
    
    core.gen_parameters(filename_parameters=filename_parameters,
                        filename=filename,
                        fileout=fileout,
                        filename_noise=filename_noise,
                        filename_init_spec=filename_init_spec,
                        filename_lsf=filename_lsf,
                        n_gauss=n_gauss,
                        lambda_amp=lambda_amp,
                        lambda_mu=lambda_mu,
                        lambda_sig=lambda_sig,
                        lambda_sig_corr_narrow=lambda_sig_corr_narrow,
                        lambda_sig_corr_broad=lambda_sig_corr_broad,
                        lambda_mu_corr_narrow=lambda_mu_corr_narrow,
                        lambda_mu_corr_broad=lambda_mu_corr_broad,
                        lambda_var_sig=lambda_var_sig,
                        lambda_r=lambda_r,
                        lb_sig=lb_sig,
                        ib_sig=ib_sig,
                        ub_sig=ub_sig,
                        lb_amp=lb_amp,
                        ub_amp=ub_amp,
                        maxiter=maxiter,
                        iprint=iprint,
                        descent=descent
                        )
    
    core.run(filename_parameters, nohup=False)
