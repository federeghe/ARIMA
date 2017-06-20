CC=gcc
CFLAGS=-c -Wall -O3
CXX=g++
CXXFLAGS=-c -Wall -O3
LDFLAGS=
CSOURCES=arima.c helpers.c
CXXSOURCES=example.cpp
CHEADERS=arima.h config.h global.h helpers.h
CXXHEADERS=
COBJECTS=$(CSOURCES:.c=.o)
CXXOBJECTS=$(CXXSOURCES:.cpp=.o)
EXECUTABLE=arimatest

all: $(CSOURCES) $(CXXSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS) $(CXXOBJECTS)
	$(CXX) $(LDFLAGS) $(COBJECTS) $(CXXOBJECTS) -o $@

.c.o: $(CHEADERS)
	$(CC) $(CFLAGS) $< -o $@

.cpp.o: $(CXXHEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(COBJECTS) $(CXXOBJECTS) $(EXECUTABLE)
