#!/bin/bash
DIR=/Users/sbastkow/Projects/Tradis/TriclosanPaper/UsedData/Triclosan/
controlDir=/Users/sbastkow/Projects/Tradis/TriclosanPaper/UsedData/BW14LB

control1=${controlDir}/controlLBrep3_S49_L999_R1_001.out.gz.CP009273.insert_site_plot.gz
control2=${controlDir}/controlLBrep4_S40_L999_R1_001.out.gz.CP009273.insert_site_plot.gz

conc=high
    echo ${conc}
    for mic in `ls ${DIR}${conc}/Rep1/combined/MIC* | awk -F'[/]' '{print $12}'`
    do
        output_dir=$(echo ${mic} | awk -F'[.]' '{print $1}')
        ./albatradis -v -a -b 0 -c 1 -q 0.01 -o ${conc}${output_dir} /Users/sbastkow/Projects/SoftwareProjects/AlbaTradis/Data/reference_BW25113.embl "${DIR}${conc}/Rep1/combined/${mic}" "${DIR}${conc}/Rep2/combined/${mic}" ${control1} ${control2}
    done;


