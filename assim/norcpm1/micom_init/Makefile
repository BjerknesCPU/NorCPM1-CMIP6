# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include Makefile.config.$(MACH)
#For bipolar grid 
CPPDEFS = -DMPI 
#for tripolar grid -DMPI -DARCTIC
#CPPDEFS = -DMPI -DDOSA 

.DEFAULT:
	-touch $@
all: micom_init
bigrid.o: bigrid.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	bigrid.F
blkdat.o: blkdat.F common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	blkdat.F
chksum.o: chksum.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	chksum.F
chksummsk.o: chksummsk.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	chksummsk.F
crc.o: crc.c crc.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c 	crc.c
dimensions.o: dimensions.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	dimensions.F
eosdat.o: eosdat.F common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	eosdat.F
eosfun.o: eosfun.F common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	eosfun.F
fill_global.o: fill_global.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	fill_global.F
geoenv_file.o: geoenv_file.F common_blocks.h common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c  geoenv_file.F
getlpi.o: getlpi.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	getlpi.F
inigeo.o: inigeo.F common_blocks.h common_geo.h mod_xc.o mod_nctools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	inigeo.F
inivar.o: inivar.F common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	inivar.F
micom_init.o: micom_init.F common_blocks.h mod_xc.o mod_nctools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c  micom_init.F
mod_nctools.o: mod_nctools.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	mod_nctools.F
mod_xc.o: mod_xc.F mod_xc_mp.inc mod_xc_sm.inc dimensions.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	mod_xc.F
pgforc.o: pgforc.F common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c 	pgforc.F
swtfrz.o: swtfrz.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	swtfrz.F
wtime.o: wtime.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	wtime.F

SRC = mod_nctools.F crc.c fill_global.F dimensions.F chksummsk.F inigeo.F eosfun.F inivar.F mod_xc.F blkdat.F geoenv_file.F pgforc.F wtime.F eosdat.F        \
swtfrz.F micom_init.F bigrid.F getlpi.F chksum.F mod_xc_mp.inc crc.h mod_xc_sm.inc common_geo.h common_blocks.h                 \
common_eos.h
OBJ = mod_nctools.o crc.o fill_global.o dimensions.o chksummsk.o inigeo.o eosfun.o inivar.o mod_xc.o blkdat.o geoenv_file.o pgforc.o wtime.o eosdat.o swtfrz.o micom_init.o bigrid.o getlpi.o chksum.o
OFF = common_geo.h bigrid.F geoenv_file.F swtfrz.F micom_init.F eosfun.F blkdat.F common_blocks.h fill_global.F crc.h          \
mod_xc_sm.inc mod_xc.F eosdat.F dimensions.F getlpi.F common_eos.h mod_xc_mp.inc inigeo.F crc.c chksummsk.F wtime.F inivar.F chksum.F mod_nctools.F          \
pgforc.F
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) ../../micom_init
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
micom_init: $(OBJ) 
	$(LD) $(OBJ) -o ../../micom_init  $(LDFLAGS)
