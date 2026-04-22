CXX = g++
CXXFLAGS = -std=c++17 -I$(HOME)/.local/include -O3 
LDFLAGS = -L$(HOME)/.local/lib -lmarqu
SRC = tfim.cpp
EXE = main

all: $(EXE) 

$(EXE): %: %.cpp
	@echo "Compiling $@ ..."
	$(CXX) $(CXXFLAGS) $< $(SRC) $(LDFLAGS) -o $@

clean:
	rm -r $(EXE) marqu_data pairs
