//------------------------------------------------------------------------------
// SPEX_CHOLDMOD/SPEX_finalize: finalize SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, 
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See 
// SPEX_CHOLMOD/License for the license. 

//------------------------------------------------------------------------------

// SPEX_finalize frees the working environment for SPEX LU library.

#include "spex_internal.h"

SPEX_info SPEX_finalize
(
    void
)
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    SPEX_mpfr_free_cache ( ) ;    // Free mpfr internal cache
    spex_gmp_finalize ( ) ;       // Reset GMP memory variables

    spex_set_initialized (false) ;
    return (SPEX_OK) ;
}

