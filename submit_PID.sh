#!/bin/bash

#copylist - Run16 Au+Au @ 200  GeV (production 2, sst+nosst stream): startLine 1, maxLine 125,618
#choose any prediefined number of lines from picoList_all_new_sst_stream.list
#coose sample form sst stream (beginning of tha file) and nosst stream (end of the file)
#sed -n '1,500 p' ./picoLists/ADD_PHYSICS_STREAM_LIST.list > picoList_submit_PID.list

#compile run macro locally and copy the compiled version (all files) in xml to scratch
starver SL17d
root -q -b -l compileRunMacroLocallyPIDeff.C

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

echo executing submitPicoHFMakerPIDeff.csh f0r picoList_submit_PID.list inside $path

csh starSubmit/submitPicoHFMakerPIDeff.csh $path ./picoLists/picoList_physics.list
