CXX = g++
CXXFLAGS = -std=c++17 -I../deps/xad/src -I../deps/xad/build/src
LDFLAGS = -L../deps/xad/build/lib
LDLIBS = -lxad
TARGET = smoke_test
SRC = main.cpp
NPROC = $(shell nproc)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC) ../deps/xad/build/lib/libxad.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(SRC) $(LDLIBS)

../deps/xad/build/lib/libxad.a:
	mkdir -p ../deps/xad/build
	cd ../deps/xad/build && cmake ..
	cd ../deps/xad/build && make -j$(NPROC)

clean:
	rm -f $(TARGET)
	rm -rf ../deps/xad/build
