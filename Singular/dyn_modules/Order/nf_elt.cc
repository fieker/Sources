#include "kernel/mod2.h" // general settings/macros
//#include"kernel/febase.h"  // for Print, WerrorS
#include"Singular/ipid.h" // for SModulFunctions, leftv
#include"Singular/number2.h" // for SModulFunctions, leftv
#include"libpolys/coeffs/numbers.h" // nRegister, coeffs.h
#include "libpolys/coeffs/coeffs.h"
#include"Singular/blackbox.h" // blackbox type
#include "nforder.h"
#include "libpolys/coeffs/bigintmat.h"

extern n_coeffType n_NFord;

static ZZ = n_InitChar(n_Z, NULL);
omBin nfElt_bin = omGetSpecBin(sizeof(nf_elt_t));

typedef nf_elt_t * nfElt_number;

static void WriteRing(const coeffs r, BOOLEAN details)
{
  ((nf *)r->data)->Print();
}

static char* CoeffString(const coeffs r)
{
  return ((nf *)r->data)->String();
}

static void EltWrite(number &a, const coeffs r)
{
  if (a) {
    nf* K = (nf*)r;
    n_Write(nf_elt_den(a), K->ord);
    StringAppendS(" / ");
    n_Write(nf_elt_num(a), ZZ);
    StringAppendS("\n");
  } else {
    StringAppendS("(Null)\n");
  }
}

number EltCreateNumDenTransfer(nf *a, number Num, number Den)
{
  nfElt_number xx = (nfElt_number) omAllocBin(nfElt_bin);
  ((bigintmat*)Num)->simplifyContentDen(&Den); //remove common factors
  nf_elt_num(xx) = Num;
  nf_elt_den(xx) = Den;
  return (number) xx;
}


number EltCreateNumDen(nf *a, number Num, number Den)
{
  Num = n_Copy(Num, nf->ord);
  Den = n_Copy(Den, ZZ);
  return EltCreateNumDenTransfer(a, Num, Den);
}


static BOOLEAN order_cmp(coeffs n, n_coeffType t, void*parameter)
{
  return (t==n_NFord) && (n->data == parameter);
}

static void KillChar(coeffs r) {
  Print("KillChar %lx\n", r);
}
#ifdef LDEBUG
  BOOLEAN EltDBTest(number, const char *, const int, const coeffs)
{
    return TRUE;
}
#endif

static void SetChar(const coeffs r)
{
  Print("%s called\n", __func__);
}
  // or NULL
  // general stuff
static number EltMult(number a, number b, const coeffs r)
{
  nforder *O = ((nf*) (r->data))->ord;
  bigintmat *c = new bigintmat(nf_elt_num(a));
  O->elMult(c, nf_elt_num(b));
  number den = n_Mult(nf_elt_den(a), nf_elt_den(b), ZZ);
  return EltCreateNumDenTransfer(r, (number)c, den);
}

static number EltAddSub(number a, number b, const coeffs r, BOOLEAN add)
{
  nforder *O = ((nf*) (r->data))->ord;
  bigintmat *c = new bigintmat((bigintmat*)nf_elt_num(a));
  number den;

  if (nf_elt_den(a) != nf_elt_den(b)) {
    den = n_Lcm(nf_elt_den(a), nf_elt_den(b), ZZ);
    number c1 = n_Div(den, nf_elt_den(a), ZZ),
    number c2 = n_Div(den, nf_elt_den(b), ZZ);
    c->skalmult(c1);
    bigintmat *d = new bigintmat((bigintmat*)nf_elt_num(b));
    d->skalmult(c2);
    if (add) 
      O->elAdd(c, (bigintmat*) d);
    else
      O->elSub(c, (bigintmat*) d);
    delete d;
    n_Delete(&c1, ZZ);
    n_Delete(&c2, ZZ);
  } else {
    den = n_Copy(nf_elt_den(a), ZZ);
    if (add) 
      O->elAdd(c, (bigintmat*) nf_elt_num(b));
    else
      O->elSub(c, (bigintmat*) nf_elt_num(b));
    delete d;
  }
  return EltCreateNumDenTransfer(r, (number)c, den);
}

static number EltSub(number a, number b, const coeffs r)
{
  return EltAddSub(a, b, r, FALSE);
}
static number EltAdd(number a, number b, const coeffs r)
{
  return EltAddSub(a, b, r, TRUE);
}
static number EltDiv(number a, number b, const coeffs r)
{
  Werror("%s called\n", __func__, a, b, r);
  return NULL;
}
static number EltIntDiv(number a, number b, const coeffs r)
{
  Werror("IntDiv called on order elts", a, b, r);
  return NULL;
}
static number EltIntMod(number a, number b, const coeffs r)
{
  Werror("IntMod called on order elts", a, b, r);
  return NULL;
}
static number EltExactDiv(number a, number b, const coeffs r)
{
  Werror("%s called\n", __func__, a, b, r);
  return NULL;
}
   /// init with an integer
static number  EltInit(long i,const coeffs r)

{
  nforder * O = (nforder*) r->data;
  if (!O) return NULL; //during init, this seems to be called with O==NULL
  coeffs C = O->basecoeffs();
  bigintmat * b = new bigintmat(O->getDim(), 1, C);
  if (O->oneIsOne()) {
    basis_elt(b, 1);
    number I = n_Init(i, C);
    b->skalmult(I, C);
    n_Delete(&I, C);
    return (number) b;
  } else
    return NULL;
}

   /// init with a GMP integer
static number  EltInitMPZ(mpz_t i, const coeffs r)

{
  Werror("%s called\n", __func__);
  return NULL;
}
   /// how complicated, (0) => 0, or positive
static int EltSize(number n, const coeffs r)

{
  Werror("%s called\n", __func__);
  return NULL;
}
   /// convertion to int, 0 if impossible
static int EltInt(number &n, const coeffs r)

{
  Werror("%s called\n", __func__);
  return NULL;
}
   /// Converts a non-negative number n into a GMP number, 0 if impossible
static void EltMPZ(mpz_t result, number &n, const coeffs r)

{
  Werror("%s called\n", __func__);
}
   /// changes argument  inline: a:= -a
   /// return -a! (no copy is returned)
   /// the result should be assigned to the original argument: e.g. a = n_Neg(a,r)
static number  EltNeg(number a, const coeffs r)
   /// return -a
{
  Werror("%s called\n", __func__);
  return NULL;
}
static number  EltInvers(number a, const coeffs r)
   /// return 1/a
{
  Werror("%s called\n", __func__);
  return NULL;
}
static number  EltCopy(number a, const coeffs r)
{
  return EltCreateMat((nforder*)r->data, (bigintmat*)a);
}

static const char * EltRead(const char * s, number * a, const coeffs r)
{
//  Print("%s called with ->%s-<\n", __func__, s);
  return s;
}

static BOOLEAN EltEqual(number a,number b, const coeffs r)
{
  Print("%s called\n", __func__, a, b, r);
  return 0;
}
static BOOLEAN EltGreater(number a,number b, const coeffs r)
{
  Print("%s called\n", __func__, a, b, r);
  return 0;
}
static BOOLEAN EltIsOne(number a, const coeffs r)
{
//  Print("%s called\n", __func__, a, r);
  return 0;
}
static BOOLEAN EltIsMOne(number a, const coeffs r)
{
//  Print("%s called\n", __func__, a, r);
  return 0;
}
static BOOLEAN EltGreaterZero(number a, const coeffs r)
{
//  Print("%s called\n", __func__, a, r);
  return 1;
}
static BOOLEAN EltIsZero(number a, const coeffs r)
{
  return (a==NULL) || ((bigintmat*)a)->isZero();
}

static nMapFunc EltSetMap(const coeffs src, const coeffs dst)
{
  Print("%s called\n", __func__, src, dst);
  return NULL;
}

static void EltDelete(number * a, const coeffs r)
{
//  Print("Deleting %lx\n%s\n", *a, (((bigintmat*)(*a))->String()));
  
  delete (bigintmat*)(*a);
  *a = NULL;
}

BOOLEAN n_nfOrderInit(coeffs r,  void * parameter)
{
  assume( getCoeffType(r) == n_NFord );
  r->nCoeffIsEqual=order_cmp;
  r->cfKillChar = KillChar;
  r->cfSetChar = SetChar;
  r->cfCoeffString=CoeffString;
  r->cfCoeffWrite=WriteRing;
  r->cfWriteShort=EltWrite;
  r->cfInit = EltInit;
  r->cfMult = EltMult;
  r->cfSub = EltSub;
  r->cfAdd = EltAdd;
  r->cfDiv = EltDiv;
  r->cfExactDiv = EltExactDiv;
  r->cfInitMPZ = EltInitMPZ;
  r->cfSize = EltSize;
  r->cfInt = EltInt;
  r->cfMPZ = EltMPZ;
  r->cfInpNeg = EltNeg;
  r->cfInvers = EltInvers;
  r->cfCopy = EltCopy;
  r->data = parameter;
  
  r->cfWriteLong = EltWrite;
  r->cfRead =EltRead;
  r->cfGreater = EltGreater;
  r->cfEqual = EltEqual;
  r->cfIsZero = EltIsZero;
  r->cfIsOne = EltIsOne;
  r->cfIsMOne = EltIsMOne;
  r->cfGreaterZero = EltGreaterZero;
  r->cfDelete = EltDelete;
  r->cfSetMap = EltSetMap;
  if (parameter)
    r->nNULL = EltInit(0, r);
#ifdef LDEBUG
  r->cfDBTest = EltDBTest;
#endif
  return FALSE;
}


