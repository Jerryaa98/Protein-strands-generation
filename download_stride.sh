wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz
mkdir stride
cp stride.tar stride/
cd stride
tar -xvf stride.tar
make
dir=$(pwd)
export PATH=$PATH:$dir
source ~/.bashrc
stride -h # for testing