//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_free: wrapper for free
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Free the memory allocated by SPEX_calloc, SPEX_malloc, or SPEX_realloc.

#include "spex_internal.h"

void SPEX_free
(
    void *p         // pointer to memory space to free
)
{
    if (!spex_initialized ( )) return ;
    SuiteSparse_free (p) ;
}

