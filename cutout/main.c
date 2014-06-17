//
//  main.c
//  cutout
//
//  Created by Eduardo on 06/02/2013.
//  Copyright (c) 2013 Institute of Astronomy. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include "fitsio.h"
#include "ast.h"

double max_array(double a[], int num_elements)
{
	int i;
	double max=-1e20;
    for (i=0; i<num_elements; i++)
    {
        if (a[i]>max)
        {
            max=a[i];
        }
    }
    return(max);
}

double min_array(double a[], int num_elements)
{
	int i;
	double min=1e20;
    for (i=0; i<num_elements; i++)
    {
        if (a[i]<min)
        {
            min=a[i];
        }
    }
    return(min);
}

int main (int argc, char *argv[]) {
	long i,j,k,l;
	FILE *pfile;
	long fsize;
	char card[FLEN_CARD];
	
	/* CFITSIO */
	fitsfile *infptr, *outfptr;
	char *inheader;
	char *options;
	int status=0, nkeys, naxis;
	long naxes_in[2] = {0,0}, naxes_out[2]={0,0}, totpix, fpixel[2], lpixel[2], inc[2];
	float *inimg, *outimg, *insubimg;
	int naxes;
	
	/* AST */
	AstFitsChan *infitschan, *outfitschan, *outfitschan0;
	AstFrameSet *inwcsinfo, *outwcsinfo, *cvt, *cvt2;
	AstFrameSet *frameseta, *framesetb;
	AstWinMap *winmap;
	const char *encode;
	
	int res=0;
	int lbnd_in[2] = {1,1}, ubnd_in[2]={0, 0};
	int lbnd_out[2] = {1,1}, ubnd_out[2]={0, 0};
	int lbnd[2] = {1,1}, ubnd[2]={0, 0};
	double ina[2], inb[2], outa[2], outb[2];
	double xin[4], yin[4], xout[4], yout[4], sx, sy;
	
	if (argc < 3) {
		printf("Description:\n");
		printf("\n");
		printf("    Extracts a cutout from an image based on a suplied header \n");
		printf("\n");
		printf("Usage:\n");
		printf("\n");
		printf("    remap image.fit stamp.hdr stamp.fit \n");
		printf("\n");
		printf("Where:\n\n");
		printf("  image.fit is the input image to extract the stamp from\n");
		printf("  stamp.hdr is the header containing the WCS of the required cutout\n");
		printf("  stamp.fit is the output cutout\n");
		printf("\n");
		return(0);
	}
	
	/* Open input file for read */
	if (fits_open_file(&infptr, argv[1], READONLY, &status)) {
		fits_report_error(stderr, status);
		return(status);
	};
	
	/* Obtain all cards in input header and concatenate them */
	if( fits_hdr2str( infptr, 0, NULL, 0, &inheader, &nkeys, &status ) ) {
		fits_report_error(stderr, status);
		return(status);
	}
    
	astBegin;
	
	/* Create a FitsChan and fill it with FITS header cards. */
	infitschan = astFitsChan( NULL, NULL, "" );
	astPutCards( infitschan, inheader );
	
	/* Read astrometry file and create a FitsChan */
	outfitschan = astFitsChan( NULL, NULL, "" );
	pfile = fopen(argv[2], "r");
	while ( fgets ( card, FLEN_CARD, pfile ) != NULL ) {
		astPutFits(outfitschan, card, 0);
	}
	fclose(pfile);
	// Rewind the FitsChan
	encode = astGetC( outfitschan, "Encoding" );
	astClear(outfitschan, "Card");
	
	/* Free the memory holding the concatenated header cards. */
	free( inheader );
    
	/* Read size of output image */
	astClear(outfitschan, "Card");
	astGetFitsI(outfitschan, "NAXIS1", &naxes);
	naxes_out[0]=naxes;
	astClear(outfitschan, "Card");
	astGetFitsI(outfitschan, "NAXIS2", &naxes);
	naxes_out[1]=naxes;
    
	/* Read WCS information from the FitsChan. */
	astClear( infitschan, "Card" );
	astClear( outfitschan, "Card" );
	inwcsinfo = astRead( infitschan );
	outwcsinfo = astRead( outfitschan );
	//astShow( outfitschan );
	
	xin[0] = 1.0;
	yin[0] = 1.0;
	xin[1] = naxes_out[0];
	yin[1] = 1.0;
	xin[2] = 1.0;
	yin[2] = naxes_out[1];
	xin[3] = naxes_out[0]/2;
	yin[3] = naxes_out[1]/2;
	astTran2(outwcsinfo, 4, xin, yin, 1, xout, yout);
	astTran2(inwcsinfo, 4, xout, yout, 0, xin, yin);
	//printf("%f %f %f %f -- %f %f\n", xin[0], xin[1], xin[2], xin[3], min_array(xin, 4), max_array(xin, 4));
	//printf("%f %f %f %f -- %f %f\n", yin[0], yin[1], yin[2], yin[3], min_array(yin, 4), max_array(yin, 4));
	
	
	frameseta = astCopy(inwcsinfo);
	framesetb = astCopy(outwcsinfo);
	
	astInvert(frameseta);
	astInvert(framesetb);
    
	cvt = astConvert(frameseta, framesetb, "SKY,PIXEL,GRID");
	if ( cvt == AST__NULL ) {
	    printf("No conversion was possible.\n");
		return(1);
    }
	
	/* Read input image */
	//fits_get_img_dim(infptr, &naxis, &status);
	fits_get_img_size(infptr, 2, naxes_in, &status);
	if (status) {
		fits_report_error(stderr, status);
		return(status);
	}
	
	totpix=naxes_in[0] * naxes_in[1];
	//inimg = (float *) malloc(totpix * sizeof(float));
	inimg = (float *) astMalloc(sizeof(float)*totpix);
	
	
	//fits_read_img (infptr, TFLOAT, 1, totpix, 0, inimg, 0, &status);
	
	//
	sx = max_array(xin, 4) - min_array(xin, 4);
	sy = max_array(yin, 4) - min_array(yin, 4);
	fpixel[0]=(long) (min_array(xin, 4) - sx/2);
	fpixel[1]=(long) (min_array(yin, 4) - sy/2);
	lpixel[0]=(long) (max_array(xin, 4) + sx/2);
	lpixel[1]=(long) (max_array(yin, 4) + sy/2);
	if (fpixel[0]<1) fpixel[0]=1;
	if (fpixel[1]<1) fpixel[1]=1;
	if (fpixel[0]>naxes_in[0]) status=1;
	if (fpixel[1]>naxes_in[1]) status=1;
    if (lpixel[0]<fpixel[0]) status=1;
	if (lpixel[1]<fpixel[1]) status=1;;
	if (lpixel[0]>naxes_in[0]) lpixel[0]=naxes_in[0];
	if (lpixel[1]>naxes_in[1]) lpixel[1]=naxes_in[1];
    
	//printf("%d %d %d %d\n", fpixel[0], fpixel[1], lpixel[0], lpixel[1]);
    
	if (!status) {
        totpix=(lpixel[0]-fpixel[0]+1)*(lpixel[1]-fpixel[1]+1);
        insubimg = (float *) astMalloc(sizeof(float)*totpix);
        inc[0]=1;
        inc[1]=1;
        fits_read_subset(infptr, TFLOAT, fpixel, lpixel, inc, 0, insubimg, 0, &status);
	}
	if (status) {
		status=0;
        //fits_report_error(stderr, status);
		//return(status);
	} else {
        
        i=0;
        for (k=fpixel[1]; k<=lpixel[1]; k++) {
            for (j=fpixel[0]; j<=lpixel[0]; j++) {
                l=j+(k-1)*(naxes_in[0]-0)-1;
                //if (i>=0) { printf("%d %d %d %d \n", i, j, k, l);}
                inimg[l]=insubimg[i];
                i++;
            }
        }
	}
	//
	
	i=0;
	for (k=fpixel[1]; k<=lpixel[1]; k++) {
		for (j=fpixel[0]; j<=lpixel[0]; j++) {
			l=j+(k-1)*(naxes_in[0]-0)-1;
			//if (i>=0) { printf("%d %d %d %d \n", i, j, k, l);}
			inimg[l]=insubimg[i];
			i++;
		}
	}
	
	//
	
	ubnd_in[0]=naxes_in[0];
	ubnd_in[1]=naxes_in[1];
	
	/* Read size of output image */
	astClear(outfitschan, "Card");
	astGetFitsI(outfitschan, "NAXIS1", &naxes);
	naxes_out[0]=naxes;
	astClear(outfitschan, "Card");
	astGetFitsI(outfitschan, "NAXIS2", &naxes);
	naxes_out[1]=naxes;
	
	totpix = naxes_out[0] * naxes_out[1];
	//outimg = (float *) malloc(totpix * sizeof(float));
	outimg = (float *) astMalloc(sizeof (float) * totpix);
	
	ubnd_out[0]=naxes_out[0];
	ubnd_out[1]=naxes_out[1];
	
	ubnd[0]=naxes_out[0];
	ubnd[1]=naxes_out[1];
	
	res = astResampleF(cvt, 2, lbnd_in, ubnd_in, inimg, NULL, AST__NEAREST, NULL, NULL, 0, 0, 500, 0.0, 2, lbnd_out, ubnd_out, lbnd, ubnd, outimg, NULL);
	
	fits_create_file(&outfptr, argv[3], &status);
	fits_create_img(outfptr, -32, 2, naxes_out, &status);
	
	astSet( outfitschan, "Card=1, Encoding=%s", encode );
	astWrite(outfitschan, outwcsinfo);
	astClear( outfitschan, "Card" );
	while ( astFindFits( outfitschan, "%f", card, 1 ) ) {
		if (strncmp(card, "NAXIS",5) == 0) continue;
		if (strncmp(card, "SIMPLE",6) == 0) continue;
		if (strncmp(card, "BITPIX",6) == 0) continue;
		fits_write_record(outfptr, card, &status);
	}
	
	fits_write_img(outfptr, TFLOAT, 1, totpix, outimg, &status);
	if (status) {
		fits_report_error(stderr, status);
		return(status);
	}
	astEnd;
	fits_close_file(outfptr, &status);
	
	return(0);
}
