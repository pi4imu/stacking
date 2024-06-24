#!/bin/bash

base=cl

xmldir=/home/aleksei/simputdir/share/sixte/instruments/srg/erosita
xml=${xmldir}/erosita_1.xml

$SIXTE/bin/runsixt \
     XMLFile=${xml} \
     RA=0.000 Dec=0.000 \
     Prefix=sim_ \
     Simput=${base}.fits \
     EvtFile=evt_${base}.fits \
     Exposure=1000
