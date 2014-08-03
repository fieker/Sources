//////////////////////////////////////////
//////////////////////////////////////////
////     Einfache Ordnungs-Klasse     ////
//////////////////////////////////////////
//////////////////////////////////////////
/*
*/
#ifndef NF_HPP
#define NF_HPP
#include "nf_elt.h"

class nf
{
private:
  ////////////////////////////////////
  ////////// Membervariablen /////////
  ////////////////////////////////////
  int rc;
  nforder * ord; // the field is going to be the quotient field of ord
                 // same basis, but elements have denominators
  coeffs R; // will be Q probably.
  
  void init(); //basic initialisation
public:
  inline int ref_count_incref(){return rc++;};
  inline int ref_count_decref(){return rc--;};
  inline int ref_count(){return rc;};
  
  
  ////////////////////////////////////
  /// 0 Konstruktoren/Destruktoren ///
  ////////////////////////////////////
  nf(nforder * ord);
  ~nf();
  void Write();
  char* String();
  void Print();

  nforder * ord(){return ord;};

  void elAdd(nf_elt_t *a, nf_elt_t *b); // a<- a+b
  void elSub(nf_elt_t *a, nf_elt_t *b); // a<- a-b
  void elMult(nf_elt_t *a, nf_elt_t *b);// a<- a*b
  void elDiv(nf_elt_t *a, nf_elt_t *b); // a<- a/b
  number elTrace(nf_elt_t *a);
  number elNorm(nf_elt_t *a);
  bigintmat * elRepMat(nf_elt_t *a);
  
};

#endif
