CC = gcc
LD = gcc
CFLAGS = -Wall -Werror
LDFLAGS = 
RM = /bin/rm -f
OBJS = january.o
EXECUTABLE = january

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE)

january.o:
	$(CC) $(CFLAGS) -c january.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)
