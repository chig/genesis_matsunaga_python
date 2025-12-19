#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8
if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana .

bin_PROGRAMS = mbar_analysis
noinst_LIBRARIES = libmbar_analysis.a
# noinst_LTLIBRARIES = libmbar_analysis.la

mbar_analysis_SOURCES = \
        mbar_matrix.fpp \
        mbar_option_str.fpp \
        mbar_option.fpp \
        mbar_control.fpp \
        mbar_setup.fpp \
        mbar_analyze.fpp \
        mbar_main.fpp

EXTRA_DIST = \
	Makefile.depends

mbar_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
mbar_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(mbar_analysis_OBJECTS) *.f90 *~ *.mod

mbar_analyze.o : mbar_analyze.fpp
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS_DISABLEFARC) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS_DISABLEFARC) -c $*.f90
endif

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
endif

mbar_analysis$(EXEEXT) : $(mbar_analysis_ADDLIBANA) $(mbar_analysis_ADDLIB) \
						$(mbar_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(mbar_analysis_OBJECTS) \
		$(mbar_analysis_ADDLIBANA) $(mbar_analysis_ADDLIB) $(LDFLAGS)

$(noinst_LIBRARIES) : $(mbar_analysis_OBJECTS)
	$(AR) $(ARFLAGS) cr $@ $(mbar_analysis_OBJECTS)

libmbar_analysis.o : ;

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
