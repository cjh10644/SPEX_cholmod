//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_calloc: wrapper for calloc
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Allocate and initialize memory space for SPEX_CHOLMOD.

#include "spex_internal.h"

void *SPEX_calloc
(
    size_t nitems,      // number of items to allocate
    size_t size         // size of each item
)
{
    if (!spex_initialized ( )) return (NULL) ;

    return (SuiteSparse_calloc (nitems, size)) ;
}

