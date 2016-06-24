#!/bin/bash

# Execute as superuser to setup Linux box for GAMESS and QMCPACK

# ===================== PACKAGES ===================== #
# install essential packages
apt-get install vim cmake git valgrind build-essential 
# for QMCPACK
apt-get install libopenmpi-dev openmpi-bin libboost-dev fftw3-dev
apt-get install liblapack-dev libhdf5-dev libxml2-dev
apt-get install python-numpy python-matplotlib
# for GAMESS
apt-get install libatlas-base-dev

# system monitoring
apt-get install lm-sensors

# ===================== .BASHRC ===================== #
# back up .bashrc
RC=~/.bashrc
if [ ! -f .bashrc_bkp ]
then cp $RC .bashrc_bkp
fi

# modify .bashrc
echo "" >> $RC
echo "# custom software"
echo "SOFTWARE=/home/`whoami`/Software" >> $RC
echo "PATH=\$PATH:\$SOFTWARE/qmcpack/build/bin:\$SOFTWARE/gamess" >> $RC
echo "PATH=\$PATH:\$SOFTWARE/qmcpack/nexus/executables" >> $RC
echo "export PYTHONPATH=\$SOFTWARE/qmcpack/nexus/library" >> $RC

# ===================== SHARED MEM ===================== #
echo "" >> /etc/sysctl.conf
echo "# increase shared memory to 2GB" >> /etc/sysctl.conf
echo "kernel.shmmax=2147483648" >> /etc/sysctl.conf
echo "kernel.shmall=2147483648" >> /etc/sysctl.conf

# ===================== PERSONAL SCRIPTS ===================== #
cd ~
git clone https://github.com/Paul-St-Young/QMC.git
git config --global user.name "Paul Young"
git config --global user.email "yyang173@illinois.edu"
git config --global core.editor vi
git config --global color.ui true
git config --global push.default simple
rm -rf Templates
mv QMC Templates
