#!/bin/bash
#need to set the directories that will be changed to XXXXX
clear

#location of offline xml path to be changed
offlinePath=/home/kswiss/Workspace/workoffline/simulator/PROD_ELVESSimulator

elvesSimulatorXML="$offlinePath/ELVESSimulator.xml"
eventFileExporterXML="$offlinePath/EventFileExporter.xml"

echo "The following files will be temporarily altered:"
echo "     $elvesSimulatorXML"
echo "     $eventFileExporterXML"

fileList="$offlinePath/ScanList2"

counter=1

#initializing the line with the placer
oldlineLON=$(echo XXXXX)
oldlineLAT=$(echo YYYYY)
oldlineOUT=$(echo ZZZZZ)

#LOOP to start the computing by simply reading the file line by line
counter=1
while IFS='' read -r line || [[ -n "$line" ]]; do

    IFS=" " read varLAT varLON <<< "$line"
    varOUT="ELVES_$counter"
    sed -i -- "s/$oldlineLAT/$varLAT/g" "$elvesSimulatorXML"
    sed -i -- "s/$oldlineLON/$varLON/g" "$elvesSimulatorXML"
    sed -i -- "s/$oldlineOUT/$varOUT/g" "$eventFileExporterXML"
    
    oldlineLAT="$varLAT"
    oldlineLON="$varLON"
    oldlineOUT="$varOUT"
    echo $oldlineOUT $oldlineLAT $oldlineLON
    (./elves);
    counter=$((counter+1))
    
done < "$fileList"


#cleaning up
sed -i -- "s/$oldlineLON/XXXXX/g" "$elvesSimulatorXML"
sed -i -- "s/$oldlineLAT/YYYYY/g" "$elvesSimulatorXML"
sed -i -- "s/$oldlineOUT/ZZZZZ/g" "$eventFileExporterXML"

echo All Cleaned and Done

