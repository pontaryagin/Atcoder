.PHONY: build bundle run

BUILD_CONFIG := RelWithDebInfo
TEST_DIR := oj-test/$(shell basename $(URL))
TARGET := Project1/main.out
BIN_PATH := ./build/$(BUILD_CONFIG)/$(TARGET)

bundle:
	cd ./Project1 && oj-bundle ./main.cpp  > main.bundle.cpp

oj-s-f: bundle
	oj submit $(URL) ./Project1/main.bundle.cpp -y

oj-s: oj-t bundle
	oj submit $(URL) ./Project1/main.bundle.cpp -y

oj-t: build
	oj test -d $(TEST_DIR) -c $(BIN_PATH)

oj-d:
	oj download $(URL) -d $(TEST_DIR)

submit:
	oj submit

build:
	cmake -B build/$(BUILD_CONFIG) -G Ninja
	cmake --build build/$(BUILD_CONFIG) --config $(BUILD_CONFIG) --target $(TARGET) --verbose

run: build
	$(BIN_PATH)

clean:
	rm -rf ./build
	cp -i ./Project1/main.template.cpp ./Project1/main.cpp
