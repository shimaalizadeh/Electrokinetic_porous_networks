# Thin wrapper around CMake. This is the default way to build:
#   make                 configure + build build/ek_solver (Release)
#   make TYPE=Debug      same build dir, Debug build type
#   make debug           Debug build in build-debug/
#   make clean | rebuild | help
#

BUILD ?= build
TYPE  ?= Release
CMAKE ?= cmake

# Use Ninja for a fresh build dir when it is installed; otherwise CMake's default.
# An existing build dir keeps the generator it was configured with.
GENERATOR := $(if $(shell command -v ninja 2>/dev/null),-G Ninja,)

.PHONY: all configure build clean rebuild debug help

all: build

configure:
	@if [ -f "$(BUILD)/CMakeCache.txt" ]; then \
	  $(CMAKE) -S . -B "$(BUILD)" -DCMAKE_BUILD_TYPE=$(TYPE); \
	else \
	  $(CMAKE) -S . -B "$(BUILD)" $(GENERATOR) -DCMAKE_BUILD_TYPE=$(TYPE); \
	fi

build: configure
	$(CMAKE) --build "$(BUILD)" --parallel
	@echo "Built $(BUILD)/ek_solver ($(TYPE))"

clean:
	rm -rf "$(BUILD)"

rebuild:
	$(MAKE) clean
	$(MAKE) all

debug:
	$(MAKE) all BUILD=build-debug TYPE=Debug

help:
	@echo "Targets:"
	@echo "  all (default)  configure + build $(BUILD)/ek_solver, TYPE=$(TYPE)"
	@echo "  clean          rm -rf $(BUILD)"
	@echo "  rebuild        clean, then all"
	@echo "  debug          Debug build in build-debug/"
	@echo "Variables: BUILD ?= build, TYPE ?= Release (Release, Debug, RelWithDebInfo)"
