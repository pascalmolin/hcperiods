#!/bin/sh
DEPS=$HOME/deps
TMP=$HOME/tmp
mkdir -p $TMP

if [[ "TRAVIS_OS_NAME" == "osx" ]]
then
    brew update
    brew install gmp mpfr
elif [[ "TRAVIS_OS_NAME" == "linux" ]] 
then
    sudo apt-get update
    sudo apt-get install libgmp-dev libmpfr-dev
else
    cd $TMP
    wget http://mpir.org/mpir-2.7.0.tar.bz2
    tar -xf mpir-2.7.0.tar.bz2
    cd mpir-2.7.0
    ./configure --enable-gmpcompat --prefix=$DEPS --disable-static
    make -j4 > /dev/null 2>&1
    make install
    GMP="--with-gmp=$DEPS"
    
    cd $TMP
    wget http://www.mpfr.org/mpfr-3.1.4/mpfr-3.1.4.tar.bz2
    tar -xf mpfr-3.1.4.tar.bz2
    cd mpfr-3.1.4
    ./configure --with-gmp=$DEPS --prefix=$DEPS --disable-static
    make -j4 > /dev/null 2>&1
    make install
    MPFR="--with-mpfr=$DEPS"
fi

cd $TMP
wget https://github.com/wbhart/flint2/archive/trunk.tar.gz
tar -xf trunk.tar.gz
cd flint2-trunk
./configure $GMP $MPFR --prefix=$DEPS --disable-static
make -j4 > /dev/null 2>&1
make install

cd $TMP
wget https://github.com/fredrik-johansson/arb/archive/2.9.0.tar.gz
tar -xf 2.9.0.tar.gz
cd arb-2.9.0
./configure $GMP $MPFR --with-flint=$DEPS --prefix=$DEPS --disable-static
make -j4 > /dev/null 2>&1
make install
