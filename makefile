default: build

.PHONY: build
build:
	-[[ -d build ]] || mkdir build
	cd build; cmake ..;  make -j8

.PHONY: format
format:
	clang-format -i src/*.cpp src/*.h src/gl/*.h src/shaders/*

.PHONY: clean
clean:
	rm -rf build/

.PHONY: run
run:
	build/bin/pls -g

