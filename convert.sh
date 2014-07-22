#!/bin/bash

molecule=$1
echo $molecule
#molecule=LiH_CISD
thres=0.01
nguess=200

pwd
convert4qmc -gamessAscii $molecule.out -ci $molecule.out -threshold $thres -readInitialGuess $nguess -add3BodyJ
