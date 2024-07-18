#!/bin/bash

base=cl
$SIXTE/bin/simputfile Simput=${base}.fits \
     Src_Name=first \
     RA=0.0 \
     Dec=0.0 \
     srcFlux=2.137e-11 \
     Elow=0.1 \
     Eup=15 \
     NBins=1000 \
     logEgrid=yes \
     Emin=2 \
     Emax=10 \
     XSPECFile=${base}.xcm
