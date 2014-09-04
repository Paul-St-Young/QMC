#!/bin/bash

# ===================== PACKAGES ===================== #
# install essential packages
apt-get install vim cmake build-essential libopenmpi-dev openmpi-bin
apt-get install libboost-dev libatlas-base-dev liblapack-dev libhdf5-dev
apt-get install fftw3-dev libxml2-dev python-numpy python-matplotlib
apt-get install git valgrind

# ===================== .BASHRC ===================== #
# back up .bashrc
RC=~/.bashrc
if [ ! -f .bashrc_bkp ]
then cp $RC .bashrc_bkp
fi

# modify .bashrc
echo "" >> $RC
echo "\# custom software"
echo "SOFTWARE=/home/`whoami`/Software" >> $RC
echo "PATH=\$PATH:\$SOFTWARE/qmcpack/build/bin:\$SOFTWARE/gamess" >> $RC
echo "PATH=\$PATH:\$SOFTWARE/qmcpack/project_suite/executables" >> $RC

# ===================== SHARED MEM ===================== #
echo "" >> /etc/sysctl.conf
echo "\# increase shared memory to 2GB" >> /etc/sysctl.conf
echo "kernel.shmmax=2147483648" >> /etc/sysctl.conf
echo "kernel.shmall=2147483648" >> /etc/sysctl.conf

# ===================== PERSONAL SCRIPTS ===================== #
cd ~
git clone https://github.com/Paul-St-Young/QMC.git
git config --global user.name "Paul Young"
git config --global user.email "yyang173@illinois.edu"
git config --global core.editor vi
git config --global color.ui true
rm -rf Templates
mv QMC Templates
