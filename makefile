CC=gcc -Wall

LIBS= -lm -lfftw3 -L/home/nadia/softwares/fftw/fftw-3.3.5/lib

all: ns2d.exe 


ns2d.exe: ns2d.o rk2.o allocate_memory.o
	$(CC) -o $@ $+ $(LIBS)


clean:
	$(RM) *.o *.exe a.out

scrub: clean
	$(RM) *.dat
