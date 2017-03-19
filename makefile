program:test.cpp 
	g++ -c -std=c++11 sequitur/symbols.cpp sequitur/symbolwrapper.cpp  StopWatch.cpp 
	g++ -pg -w -o finalopCounter -IC:\boost\boost_1_59_0 -std=c++11 test.cpp symbolwrapper.o symbols.o StopWatch.o 


	





