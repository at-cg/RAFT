EXE= raft
LIBS= -lz
CPPFLAGS= -std=c++14 -O3 

all:
	$(CXX) $(CPPFLAGS) main.cpp -o $(EXE) $(LIBS) 

split_naive:
	$(CXX) $(CPPFLAGS) split_naive.cpp -o split_naive $(LIBS)

clean:
	rm -fr $(EXE)
	rm -fr split_naive

