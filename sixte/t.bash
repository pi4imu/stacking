#!/bin/bash

### Additional part to include instrument simulation

export EXPOSURE_TIME="10000.0"
### ColNr 4 (Rrel) of "cluster.dat" times 30
export R500C_DEGREE="0.5857352413627797" # "0.190842"
### ColNr 1 (x_pix) of "cluster.dat" times 30 minus 5
export RA="356.80620847" # "17.6465"
### ColNr 2 (y_pix) of "cluster.dat" times 30 minus 5
export DEC="11.7601258" # "18.0108"

echo "Exposure time is found to be <"$EXPOSURE_TIME">"

echo "Doing instrument simulation ..."

#export LOCALFILES=`pwd`
#export HEADAS=/home/aleksei/simputdir/bin
#export LD_LIBRARY_PATH=$HEADAS/lib:$LD_LIBRARY_PATH
#export PFILES="$LOCALFILES;$HEADAS/../share/simput/pfiles;$HEADAS/../share/sixte/pfiles"
#export SIMX=$HEADAS
#export CALDB=/home/aleksei/ciao-4.15/CALDB
##$HEADAS/CALDB
#export CALDBCONFIG=$CALDB/software/tools/caldb.config
#export CALDBALIAS=$CALDB/software/tools/alias_config.fits
#export HEADASNOQUERY=
#export HEADASPROMPT=/dev/null

export HEADAS111=/home/aleksei/simputdir/bin
export FTFT=/home/aleksei/heasoft-6.30.1/heatools/x86_64-pc-linux-gnu-libc2.35/bin

echo "# Region file format: DS9 version 4.0" > ds9.reg
echo "circle(${RA},${DEC},"${R500C_DEGREE}")" >> ds9.reg

#for snr in 124 128 132 136 140; do

export snr=128

    echo "Doing all 7 eROSITA cameras for $snr ..."
    ${HEADAS111}/erosim \
	prefix="erosita_" \
	PhotonList=events_pv.fits \
	RawData=events_allpv.fits \
	background=no \
	XMLFile="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_1.xml" \
	XMLFILE1="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_1.xml" \
	XMLFILE2="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_2.xml" \
	XMLFILE3="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_3.xml" \
	XMLFILE4="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_4.xml" \
	XMLFILE5="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_5.xml" \
	XMLFILE6="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_6.xml" \
	XMLFILE7="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_7.xml" \
	Simput="../data/eROSITA_30.0x30.0/Phox/phlist_$snr.fits" \
	Exposure=$EXPOSURE_TIME \
	SkipInvalids=yes \
	seed=-1 \
	clobber=yes \
	RA=$RA \
	Dec=$DEC \
	MJDREF=50814.0
	
#done

${FTFT}/ftmerge erosita_ccd1_evt.fits,erosita_ccd2_evt.fits,erosita_ccd3_evt.fits,erosita_ccd4_evt.fits,erosita_ccd5_evt.fits,erosita_ccd6_evt.fits,erosita_ccd6_evt.fits erosita_merged_evt_$snr.fits  clobber=yes

rm erosita_ccd*.fits erosita_tel*.fits

${HEADAS111}/imgev EvtFile=erosita_merged_evt_$snr.fits Image=erosita_img_$snr.fits \
        CoordinateSystem=0 Projection=TAN  CUNIT1=deg CUNIT2=deg \
        NAXIS1=500 NAXIS2=500  CRVAL1=$RA  CRVAL2=$DEC \
        CDELT1=-0.0027778 CDELT2=0.00277778  CRPIX1=250 CRPIX2=250 \
        clobber=yes
        
${HEADAS111}/makespec EvtFile="erosita_merged_evt_$snr.fits" \
        Spectrum="erosita_spec_$snr.pha" \
        EventFilter="regfilter('ds9.reg',RA,DEC)" \
        RSPPath=${HEADAS111}/../share/sixte/instruments/srg/erosita clobber=yes
        
echo "--------------------------------------------------------------------------"


    echo "Doing all 7 eROSITA cameras for $snr ..."
    ${HEADAS111}/erosim \
	prefix="erosita_" \
	PhotonList=events_pv.fits \
	RawData=events_allpv.fits \
	background=no \
	XMLFile="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_1.xml" \
	XMLFILE1="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_1.xml" \
	XMLFILE2="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_2.xml" \
	XMLFILE3="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_3.xml" \
	XMLFILE4="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_4.xml" \
	XMLFILE5="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_5.xml" \
	XMLFILE6="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_6.xml" \
	XMLFILE7="$HEADAS111/../share/sixte/instruments/srg/erosita/erosita_7.xml" \
	Simput="../data/eROSITA_30.0x30.0/Phox/AGNphlist_$snr.fits" \
	Exposure=$EXPOSURE_TIME \
	SkipInvalids=yes \
	seed=-1 \
	clobber=yes \
	RA=$RA \
	Dec=$DEC \
	MJDREF=50814.0
	
#done

${FTFT}/ftmerge erosita_ccd1_evt.fits,erosita_ccd2_evt.fits,erosita_ccd3_evt.fits,erosita_ccd4_evt.fits,erosita_ccd5_evt.fits,erosita_ccd6_evt.fits,erosita_ccd6_evt.fits erosita_merged_evt_AGN_$snr.fits  clobber=yes

rm erosita_ccd*.fits erosita_tel*.fits

${HEADAS111}/imgev EvtFile=erosita_merged_evt_AGN_$snr.fits Image=erosita_img_AGN_$snr.fits \
        CoordinateSystem=0 Projection=TAN  CUNIT1=deg CUNIT2=deg \
        NAXIS1=500 NAXIS2=500  CRVAL1=$RA  CRVAL2=$DEC \
        CDELT1=-0.0027778 CDELT2=0.00277778  CRPIX1=250 CRPIX2=250 \
        clobber=yes
        
${HEADAS111}/makespec EvtFile="erosita_merged_evt_AGN_$snr.fits" \
        Spectrum="erosita_spec_$snr.pha" \
        EventFilter="regfilter('ds9.reg',RA,DEC)" \
        RSPPath=${HEADAS111}/../share/sixte/instruments/srg/erosita clobber=yes 
        
echo "--------------------------------------------------------------------------"

${FTFT}/ftmerge erosita_merged_evt_$snr.fits,erosita_merged_evt_AGN_$snr.fits erosita_merged_evt.fits clobber=yes 

${HEADAS111}/imgev EvtFile=erosita_merged_evt.fits Image=erosita_img.fits \
    CoordinateSystem=0 Projection=TAN  CUNIT1=deg CUNIT2=deg \
    NAXIS1=500 NAXIS2=500  CRVAL1=$RA  CRVAL2=$DEC \
    CDELT1=-0.0027778 CDELT2=0.00277778  CRPIX1=250 CRPIX2=250 \
    clobber=yes

${HEADAS111}/makespec EvtFile="erosita_merged_evt.fits" \
    Spectrum="erosita_spec.pha" \
    EventFilter="regfilter('ds9.reg',RA,DEC)" \
    RSPPath=${HEADAS111}/../share/sixte/instruments/srg/erosita clobber=yes
    
    
