# k-Means++ in C/C++

## Development Setup

Blaze uses BLAS and LAPACK libraries if they are available. Install libraries:

```bash
sudo apt-get install cmake libblas-dev liblapack-dev libboost-all-dev
```

Compile and install Blaze:

```bash
wget https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.8.tar.gz
tar -xvf blaze-3.8.tar.gz
rm -rf blaze-3.8.tar.gz
cd blaze-3.8/
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/
sudo make install
cd ..
rm -rf blaze-3.8/
```
