#include <libpolys/coeffs/bigintmat.h>
#include "nforder.h"
#include <reporter/reporter.h>
#include"libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"




////////////////////////////////////
//// Konstruktoren/Destruktoren ////
////////////////////////////////////

/*________________0_______________*/ 
void nf::init() {
  rc = 1;
  // Gibt es eine Multtable, so gibt es keine Baseorder
  ord = NULL;
  coeffs = NULL;
  // Discriminante wird erst berechnet, wenn sie benötigt wird
}

nforder::nf(nforder * O) {
  init();
  ord = O;
  coeffs = nInitChar(n_Q, 0);
}

void nf::Write() {
  StringAppendS("Field of order\n");
  ord.Write();
}

char * nf::String() {
  StringSetS("");
  Write();
  return StringEndS();
}
void nf::Print() {
  char * s = String();
  PrintS(s);
  PrintS("\n");
  omFree(s);
}
void nf_delete(nf* o) {
  if (o->ref_count_decref()>0) {
    return;
  }
  delete o;
}

nf::~nf() {
  nforder_delete(ord);
}

void nf::elAdd(nf_elt *a, nf_elt *b) {
  Werror("elAdd called");
}


void nf::elSub(nf_elt *a, nf_elt *b) {
    Werror("Error in elSub");
}

void nf::elMult(nf_elt *a, nf_elt *b) {
    Werror("Error in elMult");
}


