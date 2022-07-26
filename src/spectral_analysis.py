import itertools
import numpy as np
from casatasks import (imregrid, 
                       imsmooth,
                       importfits,
                       exportfits,
                       spxfit)
import warnings
warnings.filterwarnings("ignore")
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u


########## BEAM correction functions #######
def wsrt_beam(ang_dist):
    k = 6.6e-8
    v = 3.27e+8
    Amp = (np.cos(k*v*ang_dist))**6
    return Amp

def dist(pc, x,y):
    return np.sqrt((pc[0] - x)**2 + (pc[1] - y)**2)

def wsrt_pb_correction(wsrt_fits:str = './data/wsrt_ssp_corrected_better.fits'):
    wsrthdu = fits.open(wsrt_fits,ignore_missing_end=True)

    # Construct the primary beam for the WSRT image and divide it to make the correction
    # cos^6(kvY) where k = 6.6 × 10−8 s, v  = 327 x 10^6 Hz and Y = angular distance in radians from he pointing center
    #pixel size in radians
    wsrt_pixel_size = 22.852983*(4.848e-6)

    # total pixel
    beam_data = np.ones((1024, 1024))
    pc = 545,504 # pixel location corresponding to the pointing center

    for x,y in itertools.product(np.arange(1024),np.arange(1024)):
        pix_dist = dist(pc, x,y)
        ang_dist = pix_dist*wsrt_pixel_size
        pixel_amp = wsrt_beam(ang_dist)
        beam_data[x,y] = pixel_amp

    #Now we divide the wsrt image by the above image to correct for the beam
    wsrthdu[0].data = wsrthdu[0].data/beam_data

    #Save the beam corrected wsrt image
    wsrt_fits_pb = wsrt_fits.split('.')[0] + '_pb.fits'
    ssp_hdu = fits.PrimaryHDU (wsrthdu[0].data, header = wsrthdu[0].header)
    ssp_hdu.writeto(wsrt_fits_pb)

    # #Save the beam
    # wsrt_pb = './wsrt_pb.fits'
    # pb_hdu = fits.PrimaryHDU (beam_data, header = wsrthdu[0].header)
    # pb_hdu.writeto(wsrt_pb)

    return

########## Spectral index calculation functions #######
def read_casa(
    wsrt_fitsimage: str = './data/wsrt_ssp_corrected_better_pb.fits',
    wsrt_casaoutfile:str = './data/wsrt_ssp_corrected_pb.im',
    lofar_fitsimage:str = './data/imagefulltestfacet_caltaper60flagresidual-MFS-image-pb.fits',
    lofar_casaoutfile:str = './data/lofar_6arcsec_pb.im'
    ):
    # read images into CASA format
    importfits(
        fitsimage=wsrt_fitsimage, 
        imagename=wsrt_casaoutfile,
        overwrite=True
        )
    importfits(
        fitsimage=lofar_fitsimage, 
        imagename=lofar_casaoutfile,
        overwrite=True
        )
    

def regrid_casa(
    wsrt_casaimage: str = "./data/wsrt_ssp_corrected_pb.im",
    lofar_casaimage:str = "analysis_casa_images/lofar_6arcsec_pb.im",
    wsrt_regridoutfile:str = "./data/wsrt_ssp_corrected_pb_regridded_to_lofar_6arcsec.im",
    regrid_outfits:str = "analysis_output_fits/wsrt_ssp_corrected_pb_regridded_to_lofar_6arcsec.fits"
):
    # regrid wsrt image on LOFAR image
    imregrid(
        imagename=wsrt_casaimage, 
        output=wsrt_regridoutfile, 
        template=lofar_casaimage
        )
    # write this image to a fits
    exportfits(
        imagename=wsrt_regridoutfile, 
        fitsimage=regrid_outfits, 
        overwrite=True
        )


def smooth_casa(
    wsrt_regrid_im: str = './data/wsrt_ssp_corrected_pb_regridded_to_lofar_6arcsec.im',
    wsrt_regrid_smoothout:str = './data/wsrt_sspcorrected_pb_regrid_lofar6asec_smoothed91asecres.im',
    lofar_im:str = './data/lofar_6arcsec_pb.im',
    lofar_smoothout:str = 'analysis_casa_images/lofar_6asec_pb_smoothed91asecres.im'
    ):
    # smooth both images to a common resolution 91 arcsec (identical to the LOFAR resolution)
    # define the beam
    mybeam = {'major': '91arcsec', 'minor': '91arcsec', 'pa': '0deg'}

    # smooth the wsrt re-gridded image
    imsmooth( 
        imagename=wsrt_regrid_im, 
        kernel='gaussian', 
        beam=mybeam, 
        targetres= True, 
        outfile=wsrt_regrid_smoothout
        )

    # smooth the lofar image
    imsmooth( 
        imagename=lofar_im, 
        kernel='gaussian', 
        beam=mybeam, 
        targetres= True, 
        outfile=lofar_smoothout
        )

    # write to fits
    exportfits(
        imagename=wsrt_regrid_smoothout, 
        fitsimage="./data/wsrt_sspcorrected_pb_regrid_lofar6asec_smoothed91asecres.fits", 
        overwrite=True
    )
    exportfits(
        imagename=lofar_smoothout, 
        fitsimage="./data/lofar_6asec_pb_smoothed91asecres.fits", 
        overwrite=True
    )

def wsrt_flux_error(S, noise=0.001, bowl=-0.00046, f=0.07):
    error = np.sqrt((f*S)**2 + noise**2 + bowl**2)
    return error
def lofar_flux_error(S, noise=0.0011, f=0.10):
    error = np.sqrt((f*S)**2 + noise**2)
    return error
def alpha_error(S1, S2):
    S1_error = lofar_flux_error(S1)
    S2_error = wsrt_flux_error(S2)
    a_error = (1/(np.log(144/327)))*(np.sqrt((S1_error/S1)**2 + (S2_error/S2)**2))
    return abs(a_error)

def spectral_index(image_lofar_data, image_wsrt_data):
    listflux = [image_lofar_data, image_wsrt_data]
    freq = np.array([144,327])

    flux = np.dstack(listflux)
    log_freq = np.log10(freq)
    #log_freq = log_freq[::-1]
    print(f"Image shape: {flux.shape}")
    x,y,z = flux.shape

    alpha = np.zeros((x,y))
    alpha_e = np.zeros((x,y))

    for r in range(x):
        for c in range(y):	
            fluxvalues = flux[r,c,:]
            if not any(x==np.nan for x in fluxvalues):
                try:
                    log_fluxpix = np.log10(fluxvalues)
                    fitfunc = np.polyfit(log_freq,log_fluxpix,1)
                    #print(fitfunc[1])
                    alphapix = fitfunc[0]
                    alpha[r,c] = alphapix
                    
                    # error map
                    alpha_err = alpha_error(fluxvalues[0], fluxvalues[1])
                    alpha_e[r,c] = alpha_err
                except:
                    alpha[r,c] = 0
    print(log_fluxpix)
    print(log_freq)
    return alpha, alpha_e

def spectral_index_calculation():
    wsrt_pb_correction()
    read_casa()
    regrid_casa()
    smooth_casa()
    spectral_index_main()

def spectral_index_main():
    wsrthdu = fits.open(
        './data/wsrt_sspcorrected_pb_regrid_lofar6asec_smoothed91asecres.fits',
        ignore_missing_end=True
        )
    lofarhdu = fits.open(
        './data/lofar_6asec_pb_smoothed91asecres.fits',
        ignore_missing_end=True
        )

    lofar_data = lofarhdu[0].data[0][0]
    wsrt_data = wsrthdu[0].data

    # blank lofar data below 2 sigma
    # noise in LOFAR image in mJy
    noise_lofar = 0.00105
    lofar_data[lofar_data <= 2*noise_lofar] = np.nan
    # blank wsrt data below 2 sigma
    # noise in WSRT image in mJy
    noise_wsrt = 0.001
    wsrt_data[wsrt_data <= 2*noise_wsrt] = np.nan

    # calculate the spectral index
    spix_map, spix_err_map = spectral_index(lofar_data, wsrt_data)

    #Save the spectral index image
    spix= './data/spix_map_2sigma.fits'
    spix_hdu = fits.PrimaryHDU (spix_map, header = wsrthdu[0].header)
    spix_hdu.writeto(spix)

    #Save the spectral index image
    spix_e= './data/spix_error_map_2sigma.fits'
    spix_err_hdu = fits.PrimaryHDU (spix_err_map, header = wsrthdu[0].header)
    spix_err_hdu.writeto(spix_e)

if __name__ == "__main__":
    from argparse import ArgumentParser, RawTextHelpFormatter
    # Create parser object
    parser = ArgumentParser(
        description="Plotting script",
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument("--mode", "-m", 
                        metavar="PATH",
                        type=str, 
                        help="Re-run the entire calculation - long or just create spectral index fits file with the already run calculation results - short"
                        )
    args = parser.parse_args()

    if args.mode == "long":
        spectral_index_calculation()
    elif args.mode == "short":
        spectral_index_main()