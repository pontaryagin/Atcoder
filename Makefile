.PHONY: build bundle run

bundle:
	cd ./Project1 && oj-bundle ./main.cpp  > main.bundle.cpp

TEST_DIR := oj-test/$(shell basename $(URL))

oj-s-f: bundle
	oj s $(URL) ./Project1/main.bundle.cpp -y

oj-s: oj-t bundle
	oj s $(URL) ./Project1/main.bundle.cpp -y

oj-t:
	oj t -d $(TEST_DIR)

oj-d:
	oj d $(URL) -d $(TEST_DIR)

submit:
	oj submit

build:
	cmake -B build
	cmake --build build

run: build
	./build/Project1/main.out
