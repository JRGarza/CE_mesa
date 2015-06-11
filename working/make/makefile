#CE: copied the content from stra/work_standard_makefile here and made the necessary changes

# work_standard_makefile

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

#################################################################

ifndef STAR
STAR = star
endif

#CE: here you need to add all the additional modules/objects written for CE_mesa
# and located in WORK_SRC_DIR to be compiled
# STAR_WORK_OBJS defined in makefile_header
OBJS = CE_energy.o CE_orbit.o CE_torque.o CE_before_struct_burn_mix.o CE_after_struct_burn_mix.o $(STAR_WORK_OBJS) run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA) 
	
#################################################################



# WORK_COMPILE defined in makefile_header

run.o: $(WORK_SRC_DIR)/run.f
	$(WORK_COMPILE) $<

run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f
	$(WORK_COMPILE) $<

%.o: $(STAR_JOB_DIR)/%.f
	$(WORK_COMPILE) $<

#CE: Added a rule of how to compile the additional code that we added for CE_mesa in WORK_SRC_DIR
%.o: $(WORK_SRC_DIR)/%.f
	$(WORK_COMPILE) $<


clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
	









