CC = gcc
LD = gcc
CFLAGS = -Wall -Werror
LDFLAGS = 
RM = /bin/rm -f
OBJS = readfile.o
EXECUTABLE = readfile

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE)

readfile.o:
	$(CC) $(CFLAGS) -c readfile.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)
