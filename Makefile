EXTLIBS := -lm -lc -Llib -lport_i -lfittools -lgfortran
CCPL    := g++ -O2 -Wno-write-strings
FCPL    := gfortran -O2
FOPTS   := -fno-second-underscore -w

TARGETS := bin/example bin/recons bin/reconsHybrid bin/reconsScint             \
	bin/reconsShower

.PHONY: all clean

all: $(TARGETS)

clean:
	@$(MAKE) -C port_i clean
	@rm -rf bin lib objs $(TARGETS)

bin/%: objs/%.o lib/libfittools.a
	@mkdir -p bin
	@$(CCPL) -o $@ $< lib/libfittools.a $(EXTLIBS)
	@rm -f $<

objs/%.o: src/%.cxx src/%.h objs
	@$(CCPL) -o $@ -c $<

objs/%.o: src/%.cxx objs
	@$(CCPL) -o $@ -c $<

objs/%.o: src/%.f objs
	@$(FCPL) $(FOPTS) -o $@ -c $<

lib/libfittools.a: objs/FitTools.o objs/AS153.o lib/libport_i.a
	@mkdir -p lib
	@ar rv lib/libfittools.a objs/FitTools.o objs/AS153.o
	@ranlib lib/libfittools.a

lib/libport_i.a:
	@mkdir -p lib
	@$(MAKE) -C port_i
	@cd lib && ln -s ../port_i/lib/libport_i.a libport_i.a

objs:
	@mkdir -p objs
