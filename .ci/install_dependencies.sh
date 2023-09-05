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


# usher gets x86_64 faToVcf. So get precompiled mafft also for x86. Could
# change this later if needed by compiling mafft from source. Note - usher
# also installs mafft using at-get. We stop it from doing so because
# bash process subsitution doesn't work in that version,
# which is how ushonium python wrapper code runs it
wget https://mafft.cbrc.jp/alignment/software/mafft-7.520-linux.tgz
tar xf mafft-7.520-linux.tgz
mv mafft-linux64/mafft.bat mafft-linux64/mafft
rm mafft-7.520-linux.tgz


git clone https://github.com/yatisht/usher.git
cd usher
# There's no sudo in the singularity build
sed -i 's/sudo -E //' install/installUbuntu.sh
# Prevent mafft install because bash process subsitution doesn't work in
# the version from apt get, which is how ushonium python wrapper code runs it
sed -i 's/mafft//' install/installUbuntu.sh
# use 1 thread to compile. On my 8gb build machine, using all 4 cores dies
# because it uses all the ram
sed -i 's/make -j$(nproc)/make -j2/' install/installUbuntu.sh
./install/installUbuntu.sh
