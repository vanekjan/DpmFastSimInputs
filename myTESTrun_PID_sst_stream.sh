#!/bin/bash

#copylist - Run16 Au+Au @ 200  GeV (production 2, sst+nosst stream): startLine 1, maxLine 125,618
#choose any prediefined number of lines from picoList_all_new_sst_stream.list
#coose sample form sst stream (beginning of tha file) and nosst stream (end of the file)
sed -n '1,500 p' ./picoLists/picoList_all_new_sst_stream.list > picoList_test_sst_PID.list
sed -n '100001,100500 p' ./picoLists/picoList_all_new_sst_stream.list >> picoList_test_sst_PID.list


path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )


echo executing submitPicoHFMakerPIDeff.csh f0r picoList_test_sst_PID.list inside $path

csh starSubmit/submitPicoHFMakerPIDeff.csh $path picoList_test_sst_PID.list
