<#!/bin/bash

LOGFILE=/home/daq/GDAS/GDASscript.log
date >> $LOGFILE

TODOLIST=/home/daq/GDAS/GDAStodo.list
TODONEW=/home/daq/GDAS/GDAStodo.new

touch $TODONEW

#get current date
set -- $(date +%d%t%m%t%y)
day=$1
month_digi=$(echo $2 | sed 's/^0//')
year=$3

months=( dec jan feb mar apr may jun jul aug sep oct nov dec )
month_char=${months[$month_digi]}

#get the number of the previous week
if [ $day -le 7 ]; then
  week=5
  month_char=${months[$month_digi-1]}
  month_digi=$(($month_digi - 1))
  if (( month_digi == 2 && year != 12 )) ; then
    week=4
  fi
else
  if [ $day -le 14 ]; then
    week=1
  else
    if [ $day -le 21 ]; then
      week=2
    else
      if [ $day -le 28 ]; then
        week=3
      else
        if [ $day -gt 28 ]; then
          week=4
        fi
      fi
    fi
  fi
fi

#adding new raw file from previous week to TODO list
echo "gdas1."$month_char$year".w"$week >> $TODOLIST
echo "added gdas1."$month_char$year".w"$week" to TODO list"

LOCATION="magic1"

RAWDIR="/home/daq/GDAS/GDASraw/"
EXTDIR="/home/daq/GDAS/GDAS"$LOCATION"/"

echo "Exctracting and processing GDAS data for "$LOCATION

while read line
do
  RawFile=$RAWDIR$line
  ExtFile=$EXTDIR$line"_"$LOCATION
  echo "Processing "$line" from TODO list"

  if [ -f $ExtFile ]; then
    echo "Extracted file found, not going to process it again!"
    continue
  fi

  if [ ! -f $RawFile ]; then
    #Download the file, -c resumes broken downloads or skips if existing file is identical to the one on server
    ftp='ftp://arlftp.arlhq.noaa.gov/pub/archives/gdas1/'
    wget -c $ftp$line >> /dev/null 2>&1
    if [ -f $line ]; then
      mv $line $RawFile
    fi
  fi

  if [ -f $RawFile ]; then
    #if RawFile is found, try to extract the relevant data, unless extracted file exists
    if [ ! -f $ExtFile ]; then
      /home/daq/GDAS/bin/gdas_$LOCATION $RawFile
      mv $RAWDIR/$line"_"$LOCATION $ExtFile
      rsync -uav $ExtFile mgaug@neas.uab.cat:/media/san2/astro/CCdata/data.magic.pic.es/nfs/CCdata/M2/
    else
      echo "extracted file already exists, the data is probably in the DB, skip event and remove it from TODO list"
      echo "if the data is not in the DB yet, simply remove the file "$ExtFile" and put "$RawFile" back on the TODO list"
      continue
    fi
  else
    #if RawFile was not downloaded, put it on TODO list
    echo "raw file $RawFile not found, will be added to TODO list"
    echo $line >> $TODONEW
    continue
  fi

done < $TODOLIST


#LOCATION="magic2"
#
#RAWDIR="/home/daq/GDAS/GDASraw/"
#EXTDIR="/home/daq/GDAS/GDAS"$LOCATION"/"
#
#echo "Exctracting and processing GDAS data for "$LOCATION
#
#while read line
#do
#  RawFile=$RAWDIR$line
#  ExtFile=$EXTDIR$line"_"$LOCATION
#  echo "Processing "$line" from TODO list"
#
#  if [ -f $ExtFile ]; then
#    echo "Extracted file found, not going to process it again!"
#    continue
#  fi
#
#  if [ -f $RawFile ]; then
#    #if RawFile is found, try to extract the relevant data, unless extracted file exists
#    if [ ! -f $ExtFile ]; then
#      /home/daq/GDAS/bin/gdas_$LOCATION $RawFile
#      mv $RAWDIR/$line"_"$LOCATION $ExtFile
#    else
#      echo "extracted file already exists, the data is probably in the DB, skip event and remove it from TODO list"
#      echo "if the data is not in the DB yet, simply remove the file "$ExtFile" and put "$RawFile" back on the TODO list"
#      continue
#    fi
#  else
#    echo "raw file $RawFile not found"
#  fi
#
#done < $TODOLIST

#replace the old with the new TODO list
mv $TODONEW $TODOLIST

#remove the raw file to save space
if [ -f $RawFile ]; then
  rm -f $RawFile
fi

exit

