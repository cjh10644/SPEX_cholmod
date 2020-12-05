//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_solve: find the exact solution for Ax=b
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system LD^(-1)U x = b. It
 * essnetially serves as a wrapper for all forward and backward substitution
 * routines. This function always returns the solution matrix x as a mpz_t
 * matrix with additional pending scale factor. If a user desires to have
 * each value for each entries with pending scale factor applied, simply
 * compute x->v[j]->x[i]/x->scale and convert to double or mpfr output as
 * desired.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Unitialized on input.
 *           on output, contains the exact solution of the system
 *
 * b:        Set of RHS vectors
 *
 * A:        Input matrix. Unmodified on input/output
 *
 * L:        Lower triangular matrix. Unmodified on input/output
 *
 * U:        Upper triangular matrix. Unmodified on input/output
 *
 * S:        a dense 3*n matrix of scale factor but stored as a vector. S(1,:)
 *           is the pending scale for L, S(2,:) is the pending scale for U,
 *           and S(3, :) is the pending scale for both L and U.
 *
 * sd:       array of pivots with pending scale applied
 *
 * P, P_inv & Q_inv: the permutation vectors. unmodified on input/output.
 *
 * keep_b:   indicate if caller wants to keep the vector b. If not, this
 *           function will not make a copy of b before forward_sub and
 *           backward_sub. Instead, x->v[j]->x will be directly set as
 *           b->v[j]->x and b->v[j]->x  will be reset to NULL.
 *
 * option:   command options
 */

#define SPEX_FREE_WORK                  \
    SPEX_MPQ_CLEAR (x_scale) ;

#define SPEX_FREE_ALL                   \
    SPEX_FREE_WORK                      \
    SPEX_matrix_free (&x) ;

#include "spex_internal.h"

SPEX_info SPEX_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_matrix **x_handle, // rational solution to the system
    // input:
    SPEX_matrix *b,         // right hand side vector
    int64_t *h;             // history vector
    const SPEX_matrix *A,   // Input matrix
    const SPEX_matrix *L,   // lower triangular matrix
    const SPEX_matrix *U,   // upper triangular matrix
    const mpq_t *S,         // the pending scale factor matrix
    const mpz_t *sd,        // array of scaled pivots
    int64_t *Ldiag,         // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *P,       // row permutation
    const int64_t *P_inv,   // inverse of row permutation
    const int64_t *Q_inv,   // inverse of column permutation
    const bool keep_b,      // indicate if b will be reused
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    if (!x_handle || !A || !L || !U || !S || !sd || !P || !P_inv || !Q_inv ||
        L->m != A->m || L->n != U->m ||
        U->n != A->n || A->n != A->m || A->m != b->m )
    {
        return SPEX_INCORRECT_INPUT;
    }
    *x_handle = NULL;

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    int64_t i, j, n = L->n;
    mpq_t x_scale ;
    SPEX_MPQ_SET_NULL (x_scale) ;
    SPEX_CHECK(SPEX_mpq_init(x_scale));

    SPEX_matrix *x = NULL;   // final solution
    // even though x will be dense, we initialize it as sparse, which make each
    // of its columns initialized with length 0, and we will allocate space
    // for its mpz_t vector later depend on the value keep_b 
    SPEX_CHECK(SPEX_matrix_alloc(&x, b->n, n, true));

    for (j = 0; j < b->n; j++)
    {
        // x_scale = 1
        SPEX_CHECK(SPEX_mpq_set_ui(x_scale, 1, 1));

        //----------------------------------------------------------------------
        // solve each column of b seperately
        //----------------------------------------------------------------------
        // make a copy of b->v[j] to x->v[j]
        if (keep_b)
        {
            // just need to allocate space for x, while let i still be NULL
            x->v[j]->x = spex_create_mpz_array(n);
            if (x->v[j]->x == NULL)
            {
                SPEX_FREE_WORK;
                return SPEX_OUT_OF_MEMORY;
            }
            for (i = 0; i < n; i++)
            {
                SPEX_CHECK(SPEX_mpz_set(x->v[j]->x[i], b->v[j]->x[i]));
            }
        }
        else
        {
            x->v[j]->x = b->v[j]->x;
            b->v[j]->x = NULL;
            b->v[j]->nzmax = 0;
        }
        x->v[j]->nzmax = n;

        // bx = (LD^(-1))\b, via forward substitution
        SPEX_CHECK(spex_forward_sub(x->v[j], x_scale, h, L, Ldiag, S, sd,
            P, P_inv));

        // bx = U\bx, via back substitution
        SPEX_CHECK(spex_backward_sub(x->v[j], x_scale, U, S, sd, P, Q_inv));
    }
    //--------------------------------------------------------------------------
    // update the scale for the solution.
    //--------------------------------------------------------------------------
    // set the scaling factor x->scale = b->scale / (x_scale* A->scale)
    // the real solution is obtained by x->v[j]->x[i]/x->scale
    SPEX_CHECK(SPEX_mpq_div(x->scale, b->scale, x_scale));
    SPEX_CHECK(SPEX_mpq_div(x->scale, x->scale, A->scale));

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*x_handle) = x ;
    return (SPEX_OK) ;
}
