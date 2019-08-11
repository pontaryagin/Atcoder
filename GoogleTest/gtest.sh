echo "Compiling unit tests..."
clang++ -std=c++11 -g -Wall -Wextra -o unit_test test.cpp -lgtest -pthread
echo "Running unit tests..."
./unit_test -v
result=$?
rm -r unit_test
echo "Unit tests completed : $result"
exit $result
