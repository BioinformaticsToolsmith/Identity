# if on macos use brew to install older g++
# brew install gcc@11
#
# then use this cmake command to setup the Makefile to use the older version
cmake .. -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-11 -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-11
