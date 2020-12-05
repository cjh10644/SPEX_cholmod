//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_create_default_options: set defaults
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Create and return SPEX_options pointer with default parameters
 * upon successful allocation, which are defined in spex_internal.h
 */

#include "spex_internal.h"

SPEX_options* SPEX_create_default_options ( void )
{

    if (!spex_initialized ( )) return (NULL) ;

    //--------------------------------------------------------------------------
    // allocate the option struct
    //--------------------------------------------------------------------------

    SPEX_options* option = SPEX_malloc(sizeof(SPEX_options)) ;
    if (!option)
    {
        // out of memory
        return (NULL) ;
    }

    //--------------------------------------------------------------------------
    // set defaults
    //--------------------------------------------------------------------------

    option->print_level = SPEX_DEFAULT_PRINT_LEVEL ;
    option->prec        = SPEX_DEFAULT_PRECISION ;
    option->round       = SPEX_DEFAULT_MPFR_ROUND ;

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    return option ;
}

