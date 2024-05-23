CC = gcc
LD = gcc
CFLAGS = -Wall -Werror
LDFLAGS = 
RM = /bin/rm -f
OBJS = triang.o
EXECUTABLE = triang

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE)

triang.o:
	$(CC) $(CFLAGS) -c triang.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)
