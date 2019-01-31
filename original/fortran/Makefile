FC	= f77
FFLAGS	= -g

OBJLIBS =
LIBS    = -llapack -lblas

EXPOKIT   = expokit.o

MATAID    = mataid.o

OBJS      = clock.o $(EXPOKIT) $(MATAID) $(OBJLIBS)

goal: test.x

clean:
	rm -f *.o

test.x: $(OBJS) test.o
	$(FC) -o test.x test.o $(OBJS) $(LIBS)

.f.o:;  $(FC) $(FFLAGS) -c $<
