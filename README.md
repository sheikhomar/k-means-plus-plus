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
