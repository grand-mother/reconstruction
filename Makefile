EXTLIBS := -L/usr/lib64 -Lport_i/lib -lm -lc -lgfortran -lport_i
FOPTS   := -fno-second-underscore
CCPL    := g++ -O 
FCPL    := gfortran -O

TARGETS := recons reconsHybrid reconsScint reconsShower

.PHONY: all install clean

all: $(TARGETS) install

install: $(TARGETS)
#	@cp $(TARGETS) $(THRONG_DIR)/soft/bin/.

clean:
	rm -rf objs/*.o objs/*.a $(TARGETS)

reconsScint: objs/reconsScint.o objs/libfittools.a
	$(CCPL) -o reconsScint objs/reconsScint.o objs/libfittools.a $(EXTLIBS)

reconsShower: objs/reconsShower.o objs/libfittools.a
	$(CCPL) -o reconsShower objs/reconsShower.o objs/libfittools.a $(EXTLIBS)

recons: objs/recons.o objs/libfittools.a
	$(CCPL) -o recons objs/recons.o objs/libfittools.a $(EXTLIBS)

reconsHybrid: objs/reconsHybrid.o objs/libfittools.a
	$(CCPL) -o reconsHybrid objs/reconsHybrid.o objs/libfittools.a $(EXTLIBS)

example: objs/example.o objs/libfittools.a
	$(CCPL) -o example objs/example.o objs/libfittools.a $(EXTLIBS) 

proto6a: objs/proto6a.o objs/libfittools.a
	$(CCPL) -o proto6a objs/proto6a.o objs/libfittools.a $(EXTLIBS)

objs/recons.o: src/recons.cxx
	$(CCPL) -o objs/recons.o -c src/recons.cxx

objs/reconsHybrid.o: src/reconsHybrid.cxx
	$(CCPL) -o objs/reconsHybrid.o -c src/reconsHybrid.cxx

objs/reconsShower.o: src/reconsShower.cxx
	$(CCPL) -o objs/reconsShower.o -c src/reconsShower.cxx

objs/reconsScint.o: src/reconsScint.cxx
	$(CCPL) -o objs/reconsScint.o -c src/reconsScint.cxx

objs/example.o: src/example.cxx
	$(CCPL) -o objs/example.o -c src/example.cxx

objs/proto6a.o: src/proto6a.cxx
	$(CCPL) -o objs/proto6a.o -c src/proto6a.cxx

objs/FitTools.o: src/FitTools.h src/FitTools.cxx
	$(CCPL) -o objs/FitTools.o -c src/FitTools.cxx

objs/AS153.o: src/AS153.f
	$(FCPL) -o objs/AS153.o -c src/AS153.f

objs/libfittools.a: objs/FitTools.o objs/AS153.o
	ar rv objs/libfittools.a objs/FitTools.o objs/AS153.o
	ranlib objs/libfittools.a 
