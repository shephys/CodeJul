#!bin/bash

scons
sfmath n1=10 n2=10 output=0 >null.rsf
./sfformat den=0.5 <null.rsf > matrix.rsf
sfgrey color=j scalebar=y <matrix.rsf | sfpen


