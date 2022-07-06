IDIR=./include
ODIR=./obj
SDIR=./src

CC=nvcc -arch=native -O3 -std=c++11
CFLAGS=-I$(IDIR)

_DEPS = parameters.h structs.h util.h geom.h EOS.h output.h timestep.h init.h boundary.h killwave.h orbit.h planet.h viscosity.h solver.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = util.o geom.o EOS.o output.o timestep.o init.o boundary.o killwave.o orbit.o planet.o viscosity.o solver.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: penguin

penguin: $(OBJ)
	$(CC) -o $@ $^

$(ODIR)/main.o: main.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/util.o: $(SDIR)/util/util.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/geom.o: $(SDIR)/geom/geom.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/EOS.o: $(SDIR)/EOS/EOS.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/output.o: $(SDIR)/output/output.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/timestep.o: $(SDIR)/timestep/timestep.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/init.o: $(SDIR)/init/init.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/boundary.o: $(SDIR)/boundary/boundary.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/killwave.o: $(SDIR)/killwave/killwave.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/orbit.o: $(SDIR)/orbit/orbit.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/planet.o: $(SDIR)/planet/planet.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/viscosity.o: $(SDIR)/viscosity/viscosity.cu $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

$(ODIR)/solver.o: $(SDIR)/solver/solve.cu $(SDIR)/solver/sweeps.cu $(SDIR)/solver/riemann/* $(SDIR)/solver/recon/* $(SDIR)/solver/force/* $(SDIR)/solver/advection/* $(SDIR)/solver/dust/* $(DEPS)
	$(CC) --device-c -o $@ $< $(CFLAGS)

clean:
	rm -f penguin *~ *linkinfo
	rm -f $(ODIR)/*.o $(ODIR)/*linkinfo
	rm -f $(IDIR)/*~
	rm -f $(SDIR)/util/*~
	rm -f $(SDIR)/geom/*~
	rm -f $(SDIR)/EOS/*~
	rm -f $(SDIR)/output/*~
	rm -f $(SDIR)/init/*~
	rm -f $(SDIR)/boundary/*~
	rm -f $(SDIR)/killwave/*~
	rm -f $(SDIR)/orbit/*~
	rm -f $(SDIR)/planet/*~
	rm -f $(SDIR)/viscosity/*~
	rm -f $(SDIR)/solver/*~
	rm -f $(SDIR)/solver/recon/*~
