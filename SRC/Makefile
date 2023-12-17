TARGET = Rosenzweig.x
SHELL = /bin/sh
CPP = g++
CPPFLAGS = -O3 -Wall
LIBS = -lgsl -lgslcblas

OBJS = Rosenzweig.o

$(TARGET): $(OBJS)
	$(CPP) -o $(TARGET) $(CPPFLAGS) $^ $(LIBS)

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)

Rosenzweig.o :