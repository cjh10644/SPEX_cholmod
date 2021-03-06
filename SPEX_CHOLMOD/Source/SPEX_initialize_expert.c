//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_initialize_expert: intialize SPEX_CHOLMOD memory
// functions for GMP
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// SPEX_initialize_expert initializes the working environment for SPEX_CHOLMOD
// with custom memory functions that are used for SPEX_CHOLMOD and GMP.

// The four inputs to this function are pointers to four functions with the
// same signatures as the ANSI C malloc, calloc, realloc, and free functions.
// That is:

//     #include <stdlib.h>
//     void *malloc (size_t size) ;
//     void *calloc (size_t nmemb, size_t size) ;
//     void *realloc (void *ptr, size_t size) ;
//     void free (void *ptr) ;

#include "spex_internal.h"

SPEX_info SPEX_initialize_expert
(
    void* (*MyMalloc) (size_t),             // user-defined malloc
    void* (*MyCalloc) (size_t, size_t),     // user-defined calloc
    void* (*MyRealloc) (void *, size_t),    // user-defined realloc
    void  (*MyFree) (void *)                // user-defined free
)
{

    if (spex_initialized ( )) return (SPEX_PANIC) ;

    //--------------------------------------------------------------------------
    // define the malloc/calloc/realloc/free functions 
    //--------------------------------------------------------------------------

    SuiteSparse_config.malloc_func  = MyMalloc ;
    SuiteSparse_config.calloc_func  = MyCalloc ;
    SuiteSparse_config.realloc_func = MyRealloc ;
    SuiteSparse_config.free_func    = MyFree ;

    //--------------------------------------------------------------------------
    // Set GMP memory functions
    //--------------------------------------------------------------------------

    return (SPEX_initialize ( )) ;
}

