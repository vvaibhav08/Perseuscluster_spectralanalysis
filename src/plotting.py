import numpy as np
import aplpy
import warnings
warnings.filterwarnings("ignore")
from astropy.io import fits
from astropy.wcs import WCS

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 3.5
mpl.rcParams['xtick.major.width'] = 3.5
mpl.rcParams['font.family'] = 'serif'

def read_lofar_fits(lofar_smoothed_fits):
    """
    read_lofar_fits 
    read lofar fits and remove unneccessary header info
    remove extra dimensions from the HDU

    Parameters
    ----------
    lofar_smoothed_fits : _type_
        lofar fits file

    Returns
    -------
    HDU
        lofar fits updated HDU
    """    
    # read the lofar smoothed file for adding contours to the spix plot
    lofarhdu = fits.open(
    lofar_smoothed_fits,
    ignore_missing_end=True
    )
    # remove axis3 and axis4 header info
    lofarhdu[0].header.remove('NAXIS3')
    lofarhdu[0].header.remove('NAXIS4')
    lofarhdu[0].header.remove('CRPIX3')
    lofarhdu[0].header.remove('CRVAL3')
    lofarhdu[0].header.remove('CDELT3')
    lofarhdu[0].header.remove('CUNIT3')
    lofarhdu[0].header.remove('CTYPE3')
    lofarhdu[0].header.remove('CRPIX4')
    lofarhdu[0].header.remove('CRVAL4')
    lofarhdu[0].header.remove('CDELT4')
    lofarhdu[0].header.remove('CUNIT4')
    lofarhdu[0].header.remove('CTYPE4')
    lofarhdu[0].header.remove('PC3_1')
    lofarhdu[0].header.remove('PC4_1')
    lofarhdu[0].header.remove('PC3_2')
    lofarhdu[0].header.remove('PC4_2')
    lofarhdu[0].header.remove('PC3_3')
    lofarhdu[0].header.remove('PC4_3')
    lofarhdu[0].header.remove('PC3_4')
    lofarhdu[0].header.remove('PC4_4')
    lofarhdu[0].header.remove('PC1_3')
    lofarhdu[0].header.remove('PC2_3')
    lofarhdu[0].header.remove('PC1_4')
    lofarhdu[0].header.remove('PC2_4')
    lofarhdu[0].header['NAXIS'] = 2
    lofarhdu[0].data = lofarhdu[0].data[0][0]
    
    return lofarhdu

def read_spix_fits(spix_file:str):
    """
    read_spix_fits

    Parameters
    ----------
    spix_file : str
        path to the spectral index / spectral index error map fits file

    Returns
    -------
    HDU
        fits HDU object
    """    
    spixhdu = fits.open(
        spix_file,
        ignore_missing_end=True
        )
    # set zeros to nan
    spixhdu[0].data[spixhdu[0].data==0.000] = np.nan
    return spixhdu

def plot(spixhdu, lofarhdu, output_png:str="./spectralindex_map.png"):
    """
    plot 
    plot the given spectral index map and contours

    Parameters
    ----------
    spixhdu : HDU
        HDU for the spectral index fits
    lofarhdu : HDU
        HDU for the lofar smoothed image
    output_png : str, optional
        output png file path, by default "./spectralindex_map.png"
    """        
    fig2 = aplpy.FITSFigure(spixhdu[0], figsize=(12,12))
    cmap = mpl.cm.jet
    cmap.set_bad('white')
    fig2.show_colorscale(cmap=cmap,stretch = 'linear',vmin=-3,vmid=-0.01,vmax=0.01)
    fig2.show_contour(lofarhdu[0], levels=[-0.003,0.003,0.006,0.012,0.024],smooth=5, colors = 'black', linewidths=1)
    # fig2.recenter(49.9506671, 41.5116961, width=2, height=2)
    fig2.add_scalebar(0.71429, '1 Mpc', color='black')
    fig2.scalebar.set_font(size=20, weight='semibold', \
                        stretch='normal', family='serif', \
                        style='normal', variant='normal')
    fig2.scalebar.set_linewidth(3)
    fig2.tick_labels.set_font(size=20)
    fig2.ticks.set_linewidth(1.5)
    fig2.ticks.set_length(5)
    fig2.axis_labels.set_font(size = 20)
    fig2.add_colorbar()
    fig2.colorbar.set_font(size=15, weight='medium', \
                        stretch='normal', family='serif', \
                        style='normal', variant='normal')
    # fig2.add_grid()
    # fig2.grid.set_color('grey')
    # fig2.grid.set_alpha(0.5)
    # fig2.grid.show()
    fig2.frame.set_linewidth(2)
    fig2.save(output_png)

    return



if __name__ == "__main__":
    from argparse import ArgumentParser, RawTextHelpFormatter
    # Create parser object
    parser = ArgumentParser(
        description="Plotting script",
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument("--spectralfits", "-spix", 
                        metavar="PATH",
                        type=str, 
                        help="path to spectral index map or spectral index error map fits file"
                        )
    parser.add_argument("--lofarfits", "--lofar", 
                        metavar="PATH", 
                        type=str,  
                        help="path to lofar map fits file"
                        )
    parser.add_argument("-o", "--outputpath", 
                        metavar="PATH", 
                        type=str, 
                        default="./output.png"
                        )
    args = parser.parse_args()


    ## plot

    # lofar_smoothed_fits = './data/lofar_6asec_pb_smoothed91asecres.fits'
    # spix_file= './data/spix_map_2sigma.fits'
    lofarhdu = read_lofar_fits(args.lofarfits)
    spixhdu = read_spix_fits(args.spectralfits)
    plot(spixhdu, lofarhdu, args.outputpath)