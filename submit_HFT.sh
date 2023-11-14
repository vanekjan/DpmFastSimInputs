#!/bin/bash

#copylist - Run16 Au+Au @ 200  GeV: startLine 1, maxLine 583172 in picoList_all_new.list
#for VPDMB triggers start at line 4783
#good HFT runs start at line 166378 (RunId > 17062047)
#choose any prediefined number of lines from picoList_all_new.list
#sed -n '166378,200000 p' ./picoLists/picoList_all_new.list > picoList_submit.list

#compile run macro locally and copy the compiled version (all files) in xml to scratch
starver SL17d
root -q -b -l compileRunMacroLocallyHFT.C

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

echo executing submitPicoHFMakerHFT.csh f0r picoList_submit.list inside $path

csh starSubmit/submitPicoHFMakerHFT.csh $path ./picoLists/picoList_physics.list
