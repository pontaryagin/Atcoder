.PHONY: build bundle run test gtest

BUILD_CONFIG := RelWithDebInfo
TEST_DIR := oj-test/$(shell test ! -z $(URL) && basename $(URL))
TARGET := src/main.out
BUILD_DIR := build/$(BUILD_CONFIG)
BIN_PATH := ./$(BUILD_DIR)/$(TARGET)

bundle:
	cd ./src && oj-bundle ./main.cpp  > main.bundle.cpp

oj-s-f: bundle
	oj submit $(URL) ./src/main.bundle.cpp -y

oj-s: oj-t bundle
	oj submit $(URL) ./src/main.bundle.cpp -y

oj-t: build
	oj test -d $(TEST_DIR) -c $(BIN_PATH)

oj-d:
	oj download $(URL) -d $(TEST_DIR)

submit:
	oj submit

configure:
	cmake -B $(BUILD_DIR) -G Ninja

build:
	make configure
	cmake --build $(BUILD_DIR) --config $(BUILD_CONFIG) --target $(TARGET) --verbose

run: build
	$(BIN_PATH)

oj-verify:
	oj-verify run

gtest:
	make configure
	cmake --build $(BUILD_DIR) --config $(BUILD_CONFIG) --verbose
	./$(BUILD_DIR)/GoogleTest/unit_test

test: oj-verify
clean:
	rm -rf ./build
	cp -i ./src/main.template.cpp ./src/main.cpp
