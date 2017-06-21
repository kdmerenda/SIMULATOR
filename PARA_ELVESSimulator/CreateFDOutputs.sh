#!/bin/bash
#need to set the directories that will be changed to XXXXX
clear

#location of simulated files
directoryPath=/media/kswiss/ExtraDrive1/ELVESSimulationData/Auger3D_7km_1000m_2000pt_500ns_SIBC_0.5c_MTLE

#location of offline xml path to be changed
offlinePath=/home/kswiss/Workspace/workoffline/FORK_ELVESSimulation

elvesSimulatorXML="$offlinePath/ELVESSimulator.xml"
eventFileExporterXML="$offlinePath/EventFileExporter.xml"

echo "The following files will be temporarily altered:"
echo "     $elvesSimulatorXML"
echo "     $eventFileExporterXML"

fileList="$offlinePath/fileList"
[ -e "$fileList" ] && rm "$fileList"
find "$directoryPath"/ -name I0_1.15* | sed "s|$directoryPath/||" >> "$fileList" 
#ls "$directoryPath/I0_8*/" >> "$fileList"

#getting the amount of files to be raytraced
fileNum=$(ls -1 $directoryPath | wc -l)
echo "Number of Directories in $directoryPath: $fileNum"
counter=1

#initializing the line with the placer
oldline=$(echo XXXXX)

#LOOP to start the computing by simply reading the file line by line
while IFS='' read -r line || [[ -n "$line" ]]; do

#    echo "OldLine: $oldline, NewLine: $line"
    sed -i -- "s/$oldline/$line/g" "$elvesSimulatorXML"
    sed -i -- "s/$oldline/$line/g" "$eventFileExporterXML"
    oldline="$line"

    #run the script for the transformation with the correct path set now. 
    #the files overight each other when ran at the same time as there is no recompilation...
    #    if [ $(echo $(($counter % 2)) ) -eq 0 ]
#    then
#	echo "$line"
	(./elves);
#    else
#	echo "$line" background
#	(./elves &);
#    fi
    let "counter = $counter + 1"
    
done < "$fileList"


#cleaning up
sed -i -- "s/$oldline/XXXXX/g" "$elvesSimulatorXML"
sed -i -- "s/$oldline/XXXXX/g" "$eventFileExporterXML"
rm fileList
echo All Cleaned and Done
