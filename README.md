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

Download Enron data:

```bash
mkdir -p data/raw
curl https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/docword.enron.txt.gz \
  --output data/raw/docword.enron.txt.gz
```

Upgrade CMake:

```bash
sudo apt install build-essential libssl-dev
wget https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2.tar.gz
tar -zxvf cmake-3.20.2.tar.gz
cd cmake-3.20.2
./bootstrap
make 
sudo make install
# sudo apt-get remove cmake cmake-data
# echo "export PATH=/usr/local/share/cmake-3.20:\$PATH" >> ~/.bashrc
```

Install Boost 1.76.0:

```
sudo apt-get update
sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev
wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz 
tar -xf boost_1_76_0.tar.gz
cd boost_1_76_0
./bootstrap.sh
./b2
sudo ./b2 install
```

Build via CMake:

```bash
cmake -S standalone -B build/standalone
cmake --build build/standalone
./build/standalone/kmeans --help
```

## Downloading Datasets

Create data directories:

```bash
mkdir -p data/raw
mkdir -p data/results
```

- US Census Data (1990)

    ```bash
    curl https://archive.ics.uci.edu/ml/machine-learning-databases/census1990-mld/USCensus1990.data.txt \
        --output data/raw/USCensus1990.data.txt
    ```

- Covertype

    ```bash
    curl https://archive.ics.uci.edu/ml/machine-learning-databases/covtype/covtype.data.gz \
        --output data/raw/covtype.data.gz
    ```

- Bag of Words Datasets

    ```bash
    curl https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/docword.enron.txt.gz \
        --output data/raw/docword.enron.txt.gz
    ```
- Tower dataset

    ```bash
    curl http://homepages.uni-paderborn.de/frahling/instances/Tower.txt \
        --output data/raw/Tower.txt
    ```

## Run experiment

Install PyEnv & Python 3.8:

```bash
git clone https://github.com/pyenv/pyenv.git ~/.pyenv
cd ~/.pyenv && src/configure && make -C src
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init --path)"' >> ~/.bashrc
source ~/.bashrc
pyenv install 3.8.11
```
