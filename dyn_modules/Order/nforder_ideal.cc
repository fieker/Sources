#include <libpolys/coeffs/bigintmat.h>
#include "nforder_ideal.h"
#include "nforder.h"
#include <reporter/reporter.h>
#include "libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"




////////////////////////////////////
//// Konstruktoren/Destruktoren ////
////////////////////////////////////

/*________________0_______________*/ 
void nforder_ideal::init() {
  memset(this, 0, sizeof(this));
}

nforder_ideal::nforder_ideal(bigintmat * _basis, const coeffs O) {
  init();
  ord = O;
  basis = new bigintmat(_basis);
}

nforder_ideal::nforder_ideal(nforder_ideal *I, int) {
  init();
  ord = I->ord;
  coeffs C = ((nforder *)ord->data)->basecoeffs();
  basis = new bigintmat(I->basis);
  if (I->norm) {
    norm = n_Copy(I->norm, C);
    norm_den = n_Copy(I->norm_den, C);
  }
  if (I->min) {
    min = n_Copy(I->min, C);
    min_den = n_Copy(I->min_den, C);
  }
}

void nforder_ideal::Write() {
  coeffs C = ((nforder *)ord->data)->basecoeffs();
  StringAppend("Ideal of order");
  n_CoeffWrite(ord);
  StringAppend("\nwith basis:\n");
  basis->Write();
  if (norm) {
    StringAppendS("and norm ");
    n_Write(norm, C);
    StringAppendS(" / ");
    n_Write(norm_den, C);
    StringAppendS(" ");
  }
  if (min) {
    StringAppendS("and min ");
    n_Write(min, C);
    StringAppendS(" / ");
    n_Write(min_den, C);
    StringAppendS(" ");
  }
}

char * nforder_ideal::String() {
  StringSetS("");
  Write();
  return StringEndS();
}
void nforder_ideal::Print() {
  char * s = String();
  PrintS(s);
  PrintS("\n");
  omFree(s);
}

nforder_ideal::~nforder_ideal() {
  if (basis) delete basis;
  coeffs C = ((nforder *)ord->data)->basecoeffs();
  if (norm) {
    n_Delete(&norm, C);
    n_Delete(&norm_den, C);
  }
  if (min) {
    n_Delete(&min, C);
    n_Delete(&min_den, C);
  }
}

nforder_ideal * nf_idAdd(nforder_ideal *A, nforder_ideal *B)
{
  assume(A->order() == B->order());
  nforder * O = (nforder*) A->order()->data;
  coeffs C = O->basecoeffs();
  bigintmat * r = new bigintmat(O->getDim(), 2*O->getDim(), C);
  r->concatcol(A->viewBasis(), B->viewBasis());
  r->hnf();
  bigintmat * t1 = new bigintmat(O->getDim(), O->getDim(), C),
            * t2 = new bigintmat(t1);
  r->splitcol(t1, t2);
  delete t1;
  delete r;
  nforder_ideal *D = new nforder_ideal(t2, A->order());
  if (O->oneIsOne())
    D->setMinTransfer(t2->get(1,1), n_Init(1, C));
  D->setNormTransfer(t2->det(), n_Init(1, C));
  delete t2;
  return D;
}

