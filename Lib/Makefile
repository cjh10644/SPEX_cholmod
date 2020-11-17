#-------------------------------------------------------------------------------
# SPEX_CHOLMOD/Lib/Makefile
#-------------------------------------------------------------------------------

# To just compile the library, do:
#       make
#
# To run a demo using the library
#       cd ../Demo ; make

LIBRARY = libspexcholmod
VERSION = 1.0.0
SO_VERSION = 1

default: library

include ../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX_CHOLMOD depends on SuiteSparse_config, M, GMP, and MPFR
LDLIBS += -lsuitesparseconfig -lm -lgmp -lmpfr

C = $(CC) $(CF) -I../Include -I../SuiteSparse_config

all: install

OBJ = \
    SPEX_calloc.o SPEX_free.o spex_cppu.o spex_get_nnz_pattern.o \
    SPEX_create_default_options.o spex_get_scattered_v.o \
    spex_create_mpz_array.o SPEX_gmp.o spex_create_scattered_vector.o \
    spex_cumsum.o spex_insert_new_entry.o spex_delete_mpz_array.o \
    spex_delete_scatter_vector.o spex_ipge.o spex_dppu1.o SPEX_LUU.o \
    spex_dppu2.o SPEX_malloc.o spex_expand_vector.o SPEX_realloc.o \
    spex_finalize_and_insert_vk.o spex_triangular_solve.o spex_find_next_nz.o \
    SPEX_initialize.o SPEX_initialize_expert.o


$(OBJ): ../Include/SPEX_CHOLMOD.h ../Source/spex_internal.h Makefile

%.o: ../Source/%.c ../Include/SPEX_CHOLMOD.h ../Source/spex_internal.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------
# Create the bin
#-------------------------------------------------------------------------------

library:
	$(MAKE) install INSTALL=$(SUITESPARSE)

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $(AR_TARGET)

#-------------------------------------------------------------------------------
# install (shared C version only)
#-------------------------------------------------------------------------------

# install SPEX_CHOLMOD
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CC) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/SPEX_CHOLMOD.h $(INSTALL_INCLUDE)
	#$(CP) ../Doc/SLIP_LU_UserGuide.pdf $(INSTALL_DOC)
	#$(CP) ../README.md $(INSTALL_DOC)/SLIP_LU_README.md
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SPEX_CHOLMOD.h
	#chmod 644 $(INSTALL_DOC)/SLIP_LU_UserGuide.pdf
	#chmod 644 $(INSTALL_DOC)/SLIP_LU_README.md

# uninstall SPEX_CPEX_CHOLMOD
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SPEX_CHOLMOD.h
	$(RM) $(INSTALL_DOC)/SLIP_LU_UserGuide.pdf
	$(RM) $(INSTALL_DOC)/SLIP_LU_README.md

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

