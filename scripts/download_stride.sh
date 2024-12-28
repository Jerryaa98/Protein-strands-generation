wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz
mkdir stride
cp stride.tar.gz stride/
cd stride
tar -xvf stride.tar.gz
make
dir=$(pwd)
export PATH=$PATH:$dir
source ~/.bashrc
stride -h # for testing