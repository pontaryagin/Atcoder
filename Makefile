.PHONY: build bundle run test gtest

BUILD_CONFIG := RelWithDebInfo
TEST_DIR := oj-test/$(shell [ ! -z $(URL) ] && basename $(URL))
TARGET := src/main.out
TARGET_UNIT_TEST := GoogleTest/unit_test
BUILD_DIR := build/$(BUILD_CONFIG)
BIN_PATH := ./$(BUILD_DIR)/$(TARGET)
DL_FLAGS := $(shell [ ! -z $(AOJ) ] && echo "--system") 

bundle:
	cd ./src && oj-bundle ./main.cpp > main.bundle.cpp

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

# # For Python
# oj-s-py: oj-t-py
# 	oj submit $(URL) ./src/main.bundle.cpp -y

# oj-s-py: oj-t-py
# 	oj submit $(URL) ./src/main.bundle.cpp -y

oj-d:
	oj download $(URL) -d $(TEST_DIR) $(DL_FLAGS) || true

oj-d-f:
	oj download $(URL) -d $(TEST_DIR) -f $(DL_FLAGS)

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

test: oj-verify gtest

init-main:
	cp -i ./src/main.template.cpp ./src/main.cpp

clean:
	rm -rf ./build
	make init-main
