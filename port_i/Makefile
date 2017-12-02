EXTLIBS    := -L/usr/lib -lc -lm -lgfortran
FCPL       := gfortran -O
CCPL       := g++ -O
FOPTS      := -fno-second-underscore -O -Wpadded -Wpacked

all: example lib/libport_i.a

clean:
	rm -rf objs/*.o lib/*.a example

example: objs/example.o objs/port_i.o objs/port3f.o
	$(CCPL) -o example objs/example.o objs/port_i.o objs/port3f.o ${EXTLIBS}

objs/example.o: src/example.cxx
	$(CCPL) -o objs/example.o -c src/example.cxx

objs/port_i.o: src/port_i.h src/port_i.cxx
	$(CCPL) -o objs/port_i.o -c src/port_i.cxx

objs/port3f.o: port/port3f.f
	rm -f port/port3f.f
	cat port/*.f > port/port3f.f
	$(FCPL) $(FOPTS) -o objs/port3f.o -c port/port3f.f

lib/libport_i.a: objs/port_i.o objs/port3f.o
	ar rv lib/libport_i.a objs/port_i.o objs/port3f.o
	ranlib lib/libport_i.a
