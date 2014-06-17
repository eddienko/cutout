import _cutout

def cutout(fitsfile, hdrfile, outfile):
    res = _cutout.cutout(fitsfile, hdrfile, outfile)
    return res