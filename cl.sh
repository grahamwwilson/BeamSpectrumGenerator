#!/bin/sh

fn=test

gfortran -g -o -fPIC -static -o ${fn} ${fn}.f `cernlib -safe pawlib mathlib kernlib packlib`

exit