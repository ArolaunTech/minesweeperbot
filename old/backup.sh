rm -r old
rsync -av --exclude='old' --exclude='.git' --exclude='CMakeCache.txt' --exclude='CMakeFiles' --exclude='cmake_install.cmake' --exclude='sweepersolver' --exclude='sweepersolverold' . old