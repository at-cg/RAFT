EXE= raft
LIBS= -lz
CPPFLAGS= -std=c++14 -O3 

all:
	$(CXX) $(CPPFLAGS) main.cpp -o $(EXE) $(LIBS) 

clean:
	rm -fr $(EXE)

