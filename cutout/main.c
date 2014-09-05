//
//  main.c
//  cutout
//
//  Created by Eduardo on 06/02/2013.
//  Copyright (c) 2013 Institute of Astronomy. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"
#include "ast.h"
#include "utils.h"

#include "main.h"

#ifdef USE_OPENMP

#include <omp.h>

#endif

static int debug=1;

int cutout(char *fitsFile, char *hdrFile, char *outFile);
int cutout_list(char *driveFile) ;
int copyheaders(fitsfile *infptr, fitsfile *outfptr);
void usage();

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

void usage()
{
    printf("Description: \n");
    printf("\n");
    printf("    Extracts a cutout from an image or list of images based on a suplied header \n");
    printf("\n");
    printf("Usage:\n");
    printf("\n");
    printf("    cutout image.fit stamp.hdr stamp.fit \n\n");
    printf("    cutout driverlist\n");
    printf("\n");
    printf("Where:\n\n");
    printf("  image.fit is the input image to extract the stamp from\n");
    printf("  stamp.hdr is the header containing the WCS of the required cutout\n");
    printf("  stamp.fit is the output cutout\n");
    printf("\n");
    printf("  driverlist is a text file containing in each line the inputs as above.\n\n");
    exit(0);
}


int main (int argc, char *argv[]) {

    int c, numberInputs;
    
    while ((c = getopt (argc, argv, "hv")) != -1)
    switch(c) {
        case 'h': // print out example header
            exit(0);
        case 'v':
            break;
        default:
            usage();
    }
    
    if (argv[optind] == NULL) usage();
    
    numberInputs = argc - optind;
    
    //printf(">>>>> %d\n", numberInputs);
    
	if ((numberInputs != 1) && (numberInputs != 3)) {
        usage();
	}
    

    if (numberInputs == 3) {
        cutout(argv[optind], argv[optind+1], argv[optind+2]);
    } else {
        cutout_list(argv[optind]);
    }
    
	return(0);
	
}

int cutout_list(char *driveFile) {
    int i,j,id,nfiles;
    FILE *pfile;
    char fitsFiles[200][80], hdrFiles[200][80], outFiles[200][80], card[256];
    
    
    if (!(pfile = fopen(driveFile, "r"))) {
        printf("ERROR: File does not exist (%s)\n", driveFile);
        return 1;
    }
    
    i=0;
    while ( fgets ( card, 256, pfile ) != NULL ) {
        sscanf(card, "%s %s %s", fitsFiles[i], hdrFiles[i], outFiles[i]);
        //printf("%s : %s : %s\n", fitsFiles[i], hdrFiles[i], outFiles[i]);
		//strcpy(fitsFiles[i], card);
        i++;
	}
	fclose(pfile);
    nfiles=i;

    
    #ifdef USE_OPENMP
    
    if (debug) {
        printf ( "  C/OpenMP version\n" );
        printf ( "\n" );
        printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
        printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );
        printf ( "\n" );
    }
    double wtime;
    wtime = omp_get_wtime();
    
    #endif
    
    #ifdef USE_OPENMP
    # pragma omp parallel for default(shared) private(id)
    #endif
        for (i=0; i<nfiles; i++) {
            id = 0;
            #ifdef USE_OPENMP
            id = omp_get_thread_num ( );
            #endif
            if (debug) printf("  Processor %d :  %s -> %s\n", id, fitsFiles[i], outFiles[i]);
            cutout(fitsFiles[i], hdrFiles[i], outFiles[i]);
        }
    #ifdef USE_OPENMP
    
    wtime = omp_get_wtime ( ) - wtime;
    
    printf ( "\n" );
    printf ( "  Elapsed wall clock time = %f\n", wtime );
    
    #endif

    return(0);
}

int cutout(char *fitsFile, char *hdrFile, char *outFile) {
    int id, ii;
	long i,j,k,l;
	FILE *pfile;
	long fsize;
	char card[FLEN_CARD];
    char *fitsFileCat;
	
	/* CFITSIO */
	fitsfile *infptr, *outfptr;
	char *inheader;
	char *options;
	int nkeys, naxis, hdutype, hdunum, status=0;
	long naxes_in[2] = {0,0}, naxes_out[2]={0,0}, totpix, fpixel[2], lpixel[2], inc[2];
	float *inimg, *outimg, *insubimg = NULL;
	int naxes;
	
	/* AST */
	AstFitsChan *infitschan, *infitschan2, *outfitschan, *outfitschan0;
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
	
	status=0;

    astBegin;
    
	/* Open input file for read */
	if (fits_open_image(&infptr, fitsFile, READONLY, &status)) {
		fits_report_error(stderr, status);
		return(status);
	};
	
    /* Initialize fitsChan */
	infitschan = astFitsChan( NULL, NULL, "" );
    
    /* Read header and add only those WCS related */
    fits_get_hdrspace(infptr, &nkeys, NULL, &status);
	for (ii = 1; ii <= nkeys; ii++)  { 
          fits_read_record(infptr, ii, card, &status); /* read keyword */
		  if (strstr(card, "CD") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "CR") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "CTYPE") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "PV") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "NAX") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "WCS") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "PC") != NULL)
			  astPutFits(infitschan, card, 0);
		  else if (strstr(card, "PR") != NULL)
			  astPutFits(infitschan, card, 0);
	  }
    
    /* Read WCS information from the FitsChan. */
    astClear( infitschan, "Card" );
    inwcsinfo = astRead( infitschan );
    

	/* Read astrometry file and create a FitsChan */
	outfitschan = astFitsChan( NULL, NULL, "" );
	pfile = fopen(hdrFile, "r");
	while ( fgets ( card, FLEN_CARD, pfile ) != NULL ) {
		astPutFits(outfitschan, card, 0);
	}
	fclose(pfile);
	// Rewind the FitsChan
	encode = astGetC( outfitschan, "Encoding" );
	astClear(outfitschan, "Card");
	outwcsinfo = astRead( outfitschan );
	
    if (outwcsinfo == AST__NULL) {
        return(-1);
    }
    
	/* Read size of output image */
	astClear(outfitschan, "Card");
	astGetFitsI(outfitschan, "NAXIS1", &naxes);
	naxes_out[0]=naxes;
	astClear(outfitschan, "Card");
	astGetFitsI(outfitschan, "NAXIS2", &naxes);
	naxes_out[1]=naxes;
    
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
	fpixel[0]=(long) (min_array(xin, 4) - sx/1);
	fpixel[1]=(long) (min_array(yin, 4) - sy/1);
	lpixel[0]=(long) (max_array(xin, 4) + sx/1);
	lpixel[1]=(long) (max_array(yin, 4) + sy/1);
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
	
	res = astResampleF(cvt, 2, lbnd_in, ubnd_in, inimg, NULL, AST__LINEAR, NULL, NULL, 0, 0, 500, AST__BAD, 2, lbnd_out, ubnd_out, lbnd, ubnd, outimg, NULL);
	
	fits_create_file(&outfptr, outFile, &status);
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
	   
    // Copy headers from input fits
    fits_get_hdu_num(infptr,  &hdunum);
    copyheaders(infptr, outfptr);
    fits_movabs_hdu(infptr, 1, &hdutype, &status);
    copyheaders(infptr, outfptr);
    fits_close_file(infptr, &status);
    
    // Copy headers from catalogue if present
    fitsFileCat = str_replace(fitsFile, ".fit", "_cat.fits");
    if (fits_open_image(&infptr, fitsFileCat, READONLY, &status)) {
        status=0;
    } else {
        
        fits_movabs_hdu(infptr, hdunum, &hdutype, &status);
        copyheaders(infptr, outfptr);
        fits_close_file(infptr, &status);
    }
    
    
    fits_update_key(outfptr, TSTRING, "PROV", fitsFile, "Originating file", &status);
    
    
	fits_close_file(outfptr, &status);
    if (status) {
		fits_report_error(stderr, status);
		return(status);
	}
    
    astFree(inimg);
    astFree(outimg);
    
    astEnd;
    
	return(status);
}

int copyheaders(fitsfile *infptr, fitsfile *outfptr) {
    
    int i, nkeys, status=0;
    float fkeyval;
    char comment[80], skeyval[80];
    
    nkeys=100;
    for (i=0; i<=nkeys; i++) {
        if (updateKeys[i].keytype == 0) break;
        if (updateKeys[i].keytype == TFLOAT) {
            fits_read_key(infptr, updateKeys[i].keytype, updateKeys[i].keyname, &fkeyval, comment, &status);
            fits_update_key(outfptr, updateKeys[i].keytype, updateKeys[i].keyname, &fkeyval, comment, &status);
            if (strncmp(updateKeys[i].keyname, updateKeys[i].keytrans, 60) != 0)
                fits_update_key(outfptr, updateKeys[i].keytype, updateKeys[i].keytrans, &fkeyval, comment, &status);
        } else if (updateKeys[i].keytype == TSTRING) {
            fits_read_key(infptr, updateKeys[i].keytype, updateKeys[i].keyname, &skeyval, comment, &status);
            fits_update_key(outfptr, updateKeys[i].keytype, updateKeys[i].keyname, skeyval, comment, &status);
            //printf("%s %s\n", updateKeys[i].keyname, skeyval);
            if (strncmp(updateKeys[i].keyname, updateKeys[i].keytrans, 60) != 0)
                fits_update_key(outfptr, updateKeys[i].keytype, updateKeys[i].keytrans, skeyval, comment, &status);
        }
        status=0;
    }
    
    
    
    return(0);
}
