/* emacs edit mode for this file is -*- C++ -*- */
/* $Id: cf_algorithm.h,v 1.4 1997-09-01 08:59:53 schmidt Exp $ */

#ifndef INCL_CF_ALGORITHM_H
#define INCL_CF_ALGORITHM_H

//{{{ docu
//
// cf_algorithm.h - declarations of higher level algorithms.
//
// This header file collects declarations of most of the
// functions in factory which implement higher level algorithms
// on canonical forms (factorization, gcd, etc.).
//
// This header file corresponds to:
// cf_chinese.cc, cf_factor.cc, cf_linsys.cc, cf_resultant.cc
//
//}}}

#include <config.h>

#include "canonicalform.h"
#include "variable.h"

/*BEGINPUBLIC*/

//{{{ declarations from cf_chinese.cc
void chineseRemainder( const CanonicalForm x1, const CanonicalForm q1, const CanonicalForm x2, const CanonicalForm q2, CanonicalForm & xnew, CanonicalForm & qnew );

void chineseRemainder( const CFArray & x, const CFArray & q, CanonicalForm & xnew, CanonicalForm & qnew );
//}}}

//{{{ declarations from cf_factor.cc
CFFList factorize ( const CanonicalForm & f, bool issqrfree = false );

CFFList factorize ( const CanonicalForm & f, const Variable & alpha );

CFFList sqrFree ( const CanonicalForm & f, bool sort = false );

bool isSqrFree ( const CanonicalForm & f );
//}}}

//{{{ declarations from cf_linsys.cc
bool linearSystemSolve( CFMatrix & M );

CanonicalForm determinant( const CFMatrix & M, int n );
//}}}

//{{{ declarations from cf_resultant.cc
CFArray subResChain ( const CanonicalForm & f, const CanonicalForm & g, Variable x );
//}}}

/*ENDPUBLIC*/

#endif /* ! INCL_CF_ALGORITHM_H */
