//
//  main.h
//  cutout
//
//  Created by Eduardo on 18/06/2014.
//  Copyright (c) 2014 IoA. All rights reserved.
//

#ifndef cutout_main_h
#define cutout_main_h

#include "fitsio.h"

struct headerKeys {
    char keyname[80];
    char keytrans[80];
    int keytype;
};

struct headerKeys updateKeys[100] = {
    {"OBJECT", "OBJECT", TSTRING},
    {"TELESCOP", "TELESCOP", TSTRING},
    {"INSTRUME", "INSTRUME", TSTRING},
    {"ESO INS FILT1 NAME", "FILTER", TSTRING},
    {"ESO TEL AIRM START", "AMSTART", TFLOAT},
    {"ESO TEL AIRM END", "AMEND", TFLOAT},
    {"MJD-OBS", "MJD-OBS", TFLOAT},
    {"DATE-OBS", "DATE-OBS", TSTRING},
    {"UTC", "UTC", TFLOAT},
    {"LST", "LST", TFLOAT},
    {"EXPTIME", "EXPTIME", TFLOAT},
    {"MAGZPT", "MAGZPT", TFLOAT},
    {"SEEING", "SEEING", TFLOAT},
    {"ELLIPTIC", "ELLIPTIC", TFLOAT},
    {"STDCRMS", "STDCRMS", TFLOAT},
    {"OBSTATUS", "OBSTATUS", TSTRING},
    {"ESOGRADE", "ESOGRADE", TSTRING},
    {"APCORPK", "APCORPK", TFLOAT},
    {"APCOR1", "APCOR1", TFLOAT},
    {"APCOR2", "APCOR2", TFLOAT},
    {"APCOR3", "APCOR3", TFLOAT},
    {"APCOR4", "APCOR4", TFLOAT},
    {"APCOR5", "APCOR5", TFLOAT},
    {"APCOR6", "APCOR6", TFLOAT},
    {"APCOR7", "APCOR7", TFLOAT},
    {"MAGZPT", "MAGZPT", TFLOAT},
    {"MAGZRR", "MAGZRR", TFLOAT},
    {"NUMZPT", "NUMZPT", TFLOAT},
    {"EXTINCT", "EXTINCT", TFLOAT},
    {"STDCRMS", "STDCRMS", TFLOAT},
    {"NUMBRMS", "NUMBRMS", TFLOAT},
    {"", "", 0}
};



#endif
