#! /bin/bash
cd $(dirname $0)
echo "Compiling unit tests..."
g++  -std=c++17  -g -Wall -Wextra -o unit_test test.cpp  -lgtest -pthread
if [  $? -gt 1  ]; then
    echo "Compile failed!"
fi
echo "Running unit tests..."
./unit_test -v
result=$?
rm -r unit_test
echo "Unit tests completed : $result"
exit $result
