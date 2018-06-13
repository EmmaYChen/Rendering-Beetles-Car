# CMake generated Testfile for 
# Source directory: /Users/chenyao/Documents/GitHub/beetles-carcar/src/ext/ptex/src/tests
# Build directory: /Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(wtest "/Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests/wtest")
add_test(rtest "/Library/Python/2.7/site-packages/cmake/data/CMake.app/Contents/bin/cmake" "-DOUT=/Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests/rtest.out" "-DDATA=/Users/chenyao/Documents/GitHub/beetles-carcar/src/ext/ptex/src/tests/rtestok.dat" "-DCMD=/Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests/rtest" "-P" "/Users/chenyao/Documents/GitHub/beetles-carcar/src/ext/ptex/src/tests/compare_test.cmake")
add_test(ftest "/Library/Python/2.7/site-packages/cmake/data/CMake.app/Contents/bin/cmake" "-DOUT=/Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests/ftest.out" "-DDATA=/Users/chenyao/Documents/GitHub/beetles-carcar/src/ext/ptex/src/tests/ftestok.dat" "-DCMD=/Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests/ftest" "-P" "/Users/chenyao/Documents/GitHub/beetles-carcar/src/ext/ptex/src/tests/compare_test.cmake")
add_test(halftest "/Users/chenyao/Documents/GitHub/beetles-carcar/build/src/ext/ptex/src/tests/halftest")
