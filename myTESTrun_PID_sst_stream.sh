#!/bin/bash

#copylist - Run16 Au+Au @ 200  GeV: startLine 1, maxLine 583172 in picoList_all_new.list
#for VPDMB triggers start at line 4783
#good HFT runs start at line 166378 (RunId > 17062047)
#choose any prediefined number of lines from picoList_all_new.list
sed -n '400101,400200 p' ./picoLists/picoList_all_new_sst_stream.list > picoList_test.list


path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )


echo executing submitPicoHFMakerPIDeff.csh f0r picoList_test.list inside $path

csh starSubmit/submitPicoHFMakerPIDeff.csh $path picoList_test.list
