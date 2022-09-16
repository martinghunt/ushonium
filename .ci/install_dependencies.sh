#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get update
apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  git \
  python3-pip \
  python3-pysam \
  python3-setuptools \
  wget


if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

git clone https://github.com/yatisht/usher.git
cd usher
# There's no sudo in the singularity build
sed -i 's/sudo -E //' install/installUbuntu.sh
# use 1 thread to compile. On my 8gb build machine, using all 4 cores dies
# because it uses all the ram
sed -i 's/make -j$(nproc)/make/' install/installUbuntu.sh
./install/installUbuntu.sh
