//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_initialize: initialize SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// SPEX_initialize initializes the working evironment for SPEX_CHOLMOD.

#include "spex_internal.h"

//------------------------------------------------------------------------------
// global variable access
//------------------------------------------------------------------------------

// a global variable, but only accessible within this file.
extern bool spex_initialize_has_been_called ;

bool spex_initialize_has_been_called = false ;

bool spex_initialized ( void )
{
    return (spex_initialize_has_been_called) ;
}

void spex_set_initialized (bool s)
{
    spex_initialize_has_been_called = s ;
}

//------------------------------------------------------------------------------
// SPEX_initialize
//------------------------------------------------------------------------------

SPEX_info SPEX_initialize ( void )
{
    if (spex_initialized ( )) return (SPEX_PANIC) ;

    mp_set_memory_functions (spex_gmp_allocate, spex_gmp_reallocate,
        spex_gmp_free) ;

    spex_set_initialized (true) ;
    return (SPEX_OK) ;
}

