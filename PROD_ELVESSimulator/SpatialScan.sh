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

fileList="$offlinePath/ScanList"

counter=1

#initializing the line with the placer
oldline=$(echo XXXXX)

#LOOP to start the computing by simply reading the file line by line
while IFS='' read -r line || [[ -n "$line" ]]; do
    
    #    echo "OldLine: $oldline, NewLine: $line"
    sed -i -- "s/$oldline/$line/g" "$elvesSimulatorXML"
    sed -i -- "s/$oldline/$line/g" "$eventFileExporterXML"
    oldline="$line"
#    echo $oldline
    (./elves);
	
done < "$fileList"


#cleaning up
sed -i -- "s/$oldline/XXXXX/g" "$elvesSimulatorXML"
sed -i -- "s/$oldline/XXXXX/g" "$eventFileExporterXML"

echo All Cleaned and Done

