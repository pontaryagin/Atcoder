.PHONY: build bundle run test gtest

BUILD_CONFIG := RelWithDebInfo
TEST_DIR := oj-test/$(shell [ ! -z $(URL) ] && basename $(URL))
TARGET := src/main.out
TARGET_UNIT_TEST := GoogleTest/unit_test
BUILD_DIR := build/$(BUILD_CONFIG)
BIN_PATH := ./$(BUILD_DIR)/$(TARGET)

bundle:
	cd ./src && oj-bundle ./main.cpp  > main.bundle.cpp

ifdef FORCE
oj-s: bundle
	make oj-t || true
	oj submit $(URL) ./src/main.bundle.cpp
else
oj-s: oj-t bundle
	oj submit $(URL) ./src/main.bundle.cpp -y
endif

oj-t: build oj-d
	oj test -d $(TEST_DIR) -c $(BIN_PATH)

oj-d:
	oj download $(URL) -d $(TEST_DIR) || true

oj-d-f:
	oj download $(URL) -d $(TEST_DIR) -f

GENERATOR := $(shell ninja --version > /dev/null 2>&1 && echo "Ninja" || echo "Unix Makefiles")

configure:
	cmake -B $(BUILD_DIR) -G "$(GENERATOR)"

build:
	make configure
	cmake --build $(BUILD_DIR) --config $(BUILD_CONFIG) --target main.out --verbose

run: build
	$(BIN_PATH)

oj-verify:
	oj-verify run

gtest:
	make configure
	cmake --build $(BUILD_DIR) --config $(BUILD_CONFIG) --target unit_test --verbose
	./$(BUILD_DIR)/GoogleTest/unit_test

test: oj-verify

init-main:
	cp -i ./src/main.template.cpp ./src/main.cpp

clean:
	rm -rf ./build
	make init-main
