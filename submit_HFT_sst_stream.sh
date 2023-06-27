#!/bin/bash

#copylist - Run16 Au+Au @ 200  GeV (production 2, sst+nosst stream): startLine 1, maxLine 125,618
#choose any prediefined number of lines from picoList_all_new_sst_stream.list
sed -n '100001,125618 p' ./picoLists/picoList_all_new_sst_stream.list > picoList_submit_sst.list

#compile run macro locally and copy the compiled version (all files) in xml to scratch
starver SL17d
root -q -b -l compileRunMacroLocally.C

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )


echo executing submitPicoHFMakerHFT.csh f0r picoList_submit_sst.list inside $path

csh starSubmit/submitPicoHFMakerHFT.csh $path picoList_submit_sst.list
