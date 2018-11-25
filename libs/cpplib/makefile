NAME=heg
CC=g++
CFLAGS=-fPIC -c
OBJS=$(NAME).o
LIB=$(NAME).so

%.o : %.cpp
	$(CC) $(CFLAGS) $< -I/usr/include/python2.7

all: $(OBJS)
	$(CC) -shared -o $(LIB) $(OBJS) -lboost_python

clean:
	rm $(OBJS) $(LIB)
