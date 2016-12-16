DEPS=$HOME/deps
TMP=$HOME/tmp
mkdir -p $TMP

cd $TMP
wget https://github.com/wbhart/flint2/archive/trunk.tar.gz
tar -xf trunk.tar.gz
cd flint2-trunk
./configure --prefix=$DEPS --disable-static
make -j4 > /dev/null 2>&1
make install

cd $TMP
wget https://github.com/fredrik-johansson/arb/archive/2.9.0.tar.gz
tar -xf 2.9.0.tar.gz
cd arb-2.9.0
./configure --with-flint=$DEPS --prefix=$DEPS --disable-static
make -j4 > /dev/null 2>&1
make install
