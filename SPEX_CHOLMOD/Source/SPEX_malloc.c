//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_malloc: wrapper for malloc
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Allocate memory space for SPEX_CHOLMOD.

#include "spex_internal.h"

void *SPEX_malloc
(
    size_t size        // size of memory space to allocate
)
{
    if (!spex_initialized ( )) return (NULL) ;
    return (SuiteSparse_malloc (1, size)) ;
}

