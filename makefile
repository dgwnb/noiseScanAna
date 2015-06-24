# -*- mode: makefile -*-
#
# Makefile containing platform dependencies different projects.
MAKEFLAGS = --no-print-directory -r -s -j2


ARCH_LOC_1 := $(wildcard $(shell root-config --prefix)/test/Makefile.arch)
ARCH_LOC_2 := $(wildcard $(shell root-config --prefix)/etc/Makefile.arch)
ARCH_LOC_3 := $(wildcard $(shell root-config --prefix)/share/doc/root/test/Makefile.arch)
ifneq ($(strip $(ARCH_LOC_1)),)
  $(info Using $(ARCH_LOC_1))
  include $(ARCH_LOC_1)
else
  ifneq ($(strip $(ARCH_LOC_2)),)
    $(info Using $(ARCH_LOC_2))
    include $(ARCH_LOC_2)
  else
    ifneq ($(strip $(ARCH_LOC_3)),)
      $(info Using $(ARCH_LOC_3))
      include $(ARCH_LOC_3)
    else
      $(error Could not find Makefile.arch!)
    endif
  endif
endif


CXXFLAGS += -Wall -Wno-overloaded-virtual -Wno-unused

INCLUDES += -I./inc -I/users/yanght/include 

VPATH	= ./src ./inc ./ws

GLIBS	+= -lTMVA -lMLP
#GLIBS 	+= -L/users/yanght/lib  -lCommonFunc -lMyParticle
GLIBS	+= -lTreePlayer -lCore -lPhysics
.PHONY:
	$(OBJS)

-include $(DEPS_Submit)
Submit:		$(OBJS_Submit)
		@echo "Linking " $@
		echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@
		@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_MakeAsimovData)
MakeAsimovData: $(OBJS_MakeAsimovData)
		@echo "Linking " $@
		echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@
		@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Asymptotic)
Asymptotic: 	$(OBJS_Asymptotic)
		@echo "Linking " $@
		echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
		@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_runAsymptoticsCLs)
bin/runAsymptoticsCLs: 	$(OBJS_runAsymptoticsCLs)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_interpolation)
bin/interpolation: 	$(OBJS_interpolation)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Pvalue)
bin/Pvalue:	 	$(OBJS_Pvalue)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Pvalue_ess)
bin/Pvalue_ess:	 	$(OBJS_Pvalue_ess)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_2Dfit)
bin/2Dfit:	 	$(OBJS_2Dfit)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_nllScan)
bin/nllScan:	 	$(OBJS_nllScan)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_contour)
bin/contour:	 	$(OBJS_contour)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_runSig)
bin/runSig:	 	$(OBJS_runSig)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_plot_two)
bin/plot_two:	 	$(OBJS_plot_two)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_imp2012)
bin/imp2012:	 	$(OBJS_imp2012)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Generate)
bin/Generate:	 	$(OBJS_Generate)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Tempalte)
bin/Tempalte:	 	$(OBJS_Tempalte)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_toyq1)
bin/toyq1:	 	$(OBJS_toyq1)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Parametrization)
bin/Parametrization:	$(OBJS_Parametrization)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Threshold)
bin/Threshold:		$(OBJS_Threshold)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

-include $(DEPS_Limit)
bin/Limit:		$(OBJS_Limit)
			@echo "Linking " $@
			echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
			@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@

bin/%	: obj/%.o  obj/CommonFunc.o obj/CommonFuncDict.o
	@echo "Linking " $@
	echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
	@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@ 


obj/%Dict.o:	%.h
		@echo "Compiling $@"
		echo rootcint -f $*Dict.cxx -c $*.h $*LinkDef.h
		rootcint -f $*Dict.cxx -c inc/$*.h inc/$*LinkDef.h
		echo $(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o $*Dict.o
		$(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o obj/$*Dict.o
		rm -f $*Dict.cxx $*Dict.h 

obj/%Dict.o:	%.hh
		@echo "Compiling $@"
		echo rootcint -f $*Dict.cxx -c $*.hh $*LinkDef.h
		rootcint -f $*Dict.cxx -c inc/$*.hh inc/$*LinkDef.h
		echo $(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o $*Dict.o
		$(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o obj/$*Dict.o
		rm -f $*Dict.cxx $*Dict.h 

obj/%.o : %.cxx
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -O2 -c $< -MD -o $@ $(INCLUDES)

obj/%.o : %.cc
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -O2 -c $< -MD -o $@ $(INCLUDES)

obj/%.o : %.C
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -O2 -c $< -MD -o $@ $(INCLUDES)

clean:
	@echo "Cleaning $<..."
	rm -fr *~ obj/*.o */*~ *_Dict.* *.a 
	@echo "Done"
