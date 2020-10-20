//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_get_nnz_pattern.c: get the row-wise or column-wise nonzero
// pattern as specified.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function finds the column-wise nonzero pattern from a
// compressed-row matrix or the row-wise nonzero pattern from a
// compressed-column matrix.

#define SPEX_FREE_WORK                      \
    SPEX_FREE(&cp);                         \
    SPEX_FREE(&rowcount);

#include "spex_internal.h"

SPEX_info spex_get_nnz_pattern    // find the specified nnz pattern 
(
    // OUTPUT:                    // should be allocated upon function call
    int64_t *Cj,                  // the col index for row-wise nnz pattern,
                                  // or the row index otherwise.
    int64_t *Cp,                  // row pointers for row-wise pattern or col
                                  // pointer otherwise
    int64_t *Cx,
    // INPUT:
    bool get_col_nnz,             // if true, get the col-wise nnz pattern
    const SPEX_matrix *L,         // the target matrix
    const SPEX_matrix *U,         // the target matrix
    const SPEX_option             // command option
)
{
    SPEX_info info;
    if (!spex_initialized ( )) return (SLIP_PANIC) ;

    if (!U->i || !U->p || !Cj || !Cp || !Cx)
    {
        return SPEX_INCORRECT_INPUT;
    }

    int64_t  n  = U->n;
    int64_t Unz = U->nz;
    int64_t *Ui = U->i;
    int64_t *Up = U->p;

    int64_t *rowcount = SPEX_calloc(n, sizeof(int64_t));
    int64_t *cp       = SPEX_calloc(n, sizeof(int64_t));
    if (!rowcount || !cp)
    {
        SPEX_FREE_WORK;
        return SPEX_OUT_OF_MEMORY;
    }

    for (int i = 0 ; i < n ; i++)
    {
        L_row_offdiag[i] = -1;
    }
    for (int64_t j = 0; j < n; j++)
    {
        for (int64_t ap = 0; ap < L->v[j]->nz; ap++)
        {
            // row index
            ci = L->v[j]->i[ap];

            if (ci == P[j])  // current entry is diagonal of L(P,:)
            {
                Ldiag[j] = ap; // get the row pointer
            }
            else            // get the last row-wise off-diagonal entries
            {
                // TODO check L->v[j]->x?
                L_row_offdiag[ci] = j; // get the column index
            }
        }
    }

    for (int ap = 0 ; ap < Unz ; ap++)
    {
        rowcount [Ui [ap]]++ ;
    }

    // compute cumulative sum of rowcount to get the row pointer
    spex_cumsum(rowcount, Cp);

    for (int i = 0 ; i < n ; i++)
    {
        cp[i] = Cp[i]; // row pointer for row i
    }
    for (int j = 0 ; j < n ; i++)
    {
        for (int ap = Up[j] ; ap < Up[j+1] ; ap++)
       {
            // TODO check U->v[j]->x?
            i = Ui[ap];
            Cj[ cp[i] ] = j ;
            Cx[ cp[i] ] = ap ;
            cp[i] ++;
        }
    }

    SPEX_FREE_WORK;
    return SPEX_OK;
}
