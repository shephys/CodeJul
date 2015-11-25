#!bin/bash

scons
sfmath n1=300 n2=300 output=0 >null.rsf
./sfformat den=0.1 <null.rsf > matrix.rsf
#sfgrey color=j scalebar=y <matrix.rsf | sfpen


