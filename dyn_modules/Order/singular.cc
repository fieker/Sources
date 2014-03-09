#include "kernel/mod2.h" // general settings/macros
#include"kernel/febase.h"  // for Print, WerrorS
#include"Singular/ipid.h" // for SModulFunctions, leftv
#include"libpolys/coeffs/numbers.h" // nRegister, coeffs.h
#include "libpolys/coeffs/coeffs.h"
#include"Singular/blackbox.h" // blackbox type
#include "nforder.h"
#include "libpolys/coeffs/bigintmat.h"

static int nforder_type_id=0;
static n_coeffType nforder_type =n_CF;

static void WriteRing(const coeffs r, BOOLEAN details)
{
  printf("%s%d:%s: %s\n", __FILE__, __LINE__, __func__, "RING");
}

static char* CoeffString(const coeffs r)
{
  return "Ring";
}
static void EltWrite(number &a, const coeffs r)
{
  bigintmat * b = (bigintmat*)a;
  if (a) {
    char * m = b->String();
    ::Print("%s\n", m);
  } else {
    Print("(Null)\n");
  }

  printf("%s%d:%s: %s\n", __FILE__, __LINE__, __func__, "RING-ELT");
}

static number EltCreateMat(nforder *a, bigintmat *b)
{
  number xx = (number) new bigintmat((bigintmat*)b);
  Print("Created new element %lx from %lx\n", xx, b);
  return (number) xx;
}


static BOOLEAN order_cmp(coeffs n, n_coeffType t, void*parameter)
{
  return (t==nforder_type) && (n->data == parameter);
}

//static void KillChar(coeffs r); //  undo all initialisations
                                // or NULL
static void SetChar(const coeffs r)
{
  Print("%s called\n", __func__);
}
                                // or NULL
   // general stuff
static number EltMult(number a, number b, const coeffs r)
{
  nforder *O = (nforder*) (r->data);
  bigintmat *c = new bigintmat((bigintmat*)a);
  O->elMult(c, (bigintmat*) b);
  return (number) c;
}
static number EltSub(number a, number b, const coeffs r)
{
  nforder *O = (nforder*) (r->data);
  bigintmat *c = new bigintmat((bigintmat*)a);
  O->elSub(c, (bigintmat*) b);
  return (number) c;
}
static number EltAdd(number a, number b, const coeffs r)
{
  nforder *O = (nforder*) (r->data);
  bigintmat *c = new bigintmat((bigintmat*)a);
  O->elAdd(c, (bigintmat*) b);
  return (number) c;
}
static number EltDiv(number a, number b, const coeffs r)
{
  Werror("%s called\n", __func__);
  return NULL;
}
static number EltIntDiv(number a, number b, const coeffs r)
{
  Werror("IntDiv called on order elts");
  return NULL;
}
static number EltIntMod(number a, number b, const coeffs r)
{
  Werror("IntMod called on order elts");
  return NULL;
}
static number EltExactDiv(number a, number b, const coeffs r)
{
  Werror("%s called\n", __func__);
  return NULL;
}
   /// init with an integer
static number  EltInit(long i,const coeffs r)

{
  Print("%s called\n", __func__);
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
   /// return 1/a
{
  Werror("%s called\n", __func__);
  return NULL;
}
static number  EltInvers(number a, const coeffs r)
   /// return a copy of a
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
  Print("%s called with ->%s-<\n", __func__, s);
  return s;
}

static BOOLEAN EltEqual(number a,number b, const coeffs r)
{
  Print("%s called\n", __func__);
  return 0;
}
static BOOLEAN EltGreater(number a,number b, const coeffs r)
{
  Print("%s called\n", __func__);
  return 0;
}
static BOOLEAN EltIsOne(number a, const coeffs r)
{
  Print("%s called\n", __func__);
  return 0;
}
static BOOLEAN EltIsMOne(number a, const coeffs r)
{
  Print("%s called\n", __func__);
  return 0;
}
static BOOLEAN EltGreaterZero(number a, const coeffs r)
{
  Print("%s called\n", __func__);
  return 0;
}
static BOOLEAN EltIsZero(number a, const coeffs r)
{
  Print("%s called\n", __func__);
  return 0;
}
static void  EltPowerSmall(number a, int i, number * result, const coeffs r)
{
  Print("%s called\n", __func__);
  result = NULL;
}

static nMapFunc EltSetMap(const coeffs src, const coeffs dst)
{
  Print("%s called\n", __func__);
  return NULL;
}

static void EltDelete(number * a, const coeffs r)
{
  Print("Deleting %lx\n%s\n", *a, (((bigintmat*)(*a))->String()));
  delete (bigintmat*)(*a);
  *a = NULL;
}

BOOLEAN n_nfOrderInit(coeffs r,  void * parameter)
{
  puts("nfOrderInit called");
  assume( getCoeffType(r) == nforder_type );
  r->nCoeffIsEqual=order_cmp;
  r->cfKillChar = NULL;
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
  r->cfNeg = EltNeg;
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
  r->cfPower = EltPowerSmall;
  r->cfDelete = EltDelete;
  r->cfSetMap = EltSetMap;
  return FALSE;
}

// coeffs stuff: -----------------------------------------------------------
static coeffs nforder_AE=NULL;
static void nforder_Register()
{
  puts("nforder_Register called");
  nforder_type=nRegister(n_unknown,n_nfOrderInit);
  nforder_AE=nInitChar(nforder_type,NULL);
}
// black box stuff: ---------------------------------------------------------
static void * nforder_Init(blackbox */*b*/)
{
  nforder_AE->ref++;
  return nforder_AE;
}
static char * nforder_String(blackbox *b, void *d)
{
  StringSetS("");
  coeffs c=(coeffs)d;
  if (c==NULL) StringAppendS("oo");
  else StringAppend("Coeff(%d)",c->type);
  Print("d: %lx\n", d);
  if (d && c->data) ((nforder *)(c->data))->Print();
  return omStrDup(StringEndS());
}
static void * nforder_Copy(blackbox*b, void *d)
{  coeffs c=(coeffs)d; if (c!=NULL) c->ref++;return d; }

static BOOLEAN nforder_Assign(leftv l, leftv r)
{
  if (l->Typ()==r->Typ())
  {
    coeffs lc=(coeffs)l->Data();
    if (lc!=NULL) lc->ref--;
    coeffs rc=(coeffs)r->Data();
    if (rc!=NULL) rc->ref++;
    if (l->rtyp==IDHDL)
    {
      IDDATA((idhdl)l->data)=(char *)rc;
    }
    else
    {
      l->data=(void *)rc;
    }
    return FALSE;
  }
  return TRUE;
}
static void nforder_destroy(blackbox *b, void *d)
{
  if (d!=NULL)
  {
    coeffs c=(coeffs)d;
    c->ref--;
  }
}
static BOOLEAN nforder_bb_setup()
{
  blackbox *b=(blackbox*)omAlloc0(sizeof(blackbox));
  // all undefined entries will be set to default in setBlackboxStuff
  // the default Print is quite usefule,
  // all other are simply error messages
  b->blackbox_destroy=nforder_destroy;
  b->blackbox_String=nforder_String;
  //b->blackbox_Print=blackbox_default_Print;
  b->blackbox_Init=nforder_Init;
  b->blackbox_Copy=nforder_Copy;
  b->blackbox_Assign=nforder_Assign;
  //b->blackbox_Op1=blackbox_default_Op1;
  //b->blackbox_Op2=blackbox_default_Op2;
  //b->blackbox_Op3=blackbox_default_Op3;
  //b->blackbox_OpM=blackbox_default_OpM;
  nforder_type_id = setBlackboxStuff(b,"NFOrder");
  Print("setup: created a blackbox type [%d] '%s'",nforder_type_id, getBlackboxName(nforder_type_id));
  PrintLn();
  return FALSE; // ok, TRUE = error!
}

// module stuff: ------------------------------------------------------------

static BOOLEAN build_ring(leftv result, leftv arg)
{

  int dimension = (int)(long)arg->Data();

  bigintmat **multtable = (bigintmat**)omAlloc(dimension*sizeof(bigintmat*));
  arg = arg->next;
  for (int i=0; i<dimension; i++) {
    multtable[i] = new bigintmat((bigintmat*)arg->Data());
    arg = arg->next;
  }
  nforder *o = new nforder(dimension, multtable, nInitChar(n_Z, 0));
  result->rtyp=nforder_type_id; // set the result type
  result->data=(char*)nInitChar(nforder_type, o);// set the result data
  Print("result is %lx\n", result->data);
  currRing->cf = (coeffs)result->data;

  return FALSE;
}

static BOOLEAN elt_from_mat(leftv result, leftv arg)
{
  coeffs C = (coeffs)arg->Data();
  nforder *O = (nforder*) (C->data);
  arg = arg->next;
  bigintmat *b = (bigintmat*) arg->Data();
  result->rtyp = NUMBER_CMD;
  result->data = (char*)EltCreateMat(O, b);
  return FALSE;
}

static BOOLEAN _calcdisc(leftv result, leftv arg)
{
  assume (arg->Typ()==nforder_type_id);
  coeffs c = (coeffs)arg->Data();
  assume (c->type = nforder_type);
  nforder * o = (nforder*)c->data;
  o->calcdisc();
  result->rtyp=NONE;

  return FALSE;
}

static BOOLEAN pMaximalOrder(leftv result, leftv arg)
{
  assume (arg->Typ()==nforder_type_id);
  coeffs c = (coeffs)arg->Data();
  assume (c->type = nforder_type);
  nforder * o = (nforder*)c->data;
  arg = arg->next;
  long p = (int)(long)arg->Data();
  number P = n_Init(p, o->basecoeffs());

  nforder *op = pmaximal(o, P);

  result->rtyp=nforder_type_id; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data
  Print("result is %lx\n", result->data);
  currRing->cf = (coeffs)result->data;

  return FALSE;
}

static BOOLEAN oneStep(leftv result, leftv arg)
{
  assume (arg->Typ()==nforder_type_id);
  coeffs c = (coeffs)arg->Data();
  assume (c->type = nforder_type);
  nforder * o = (nforder*)c->data;
  arg = arg->next;
  long p = (int)(long)arg->Data();
  number P = n_Init(p, o->basecoeffs());

  nforder *op = onestep(o, P, o->basecoeffs());

  result->rtyp=nforder_type_id; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data
  Print("result is %lx\n", result->data);
  currRing->cf = (coeffs)result->data;

  return FALSE;
}

static BOOLEAN nforder_simplify(leftv result, leftv arg)
{
  assume (arg->Typ()==nforder_type_id);
  coeffs c = (coeffs)arg->Data();
  assume (c->type = nforder_type);
  nforder * o = (nforder*)c->data;

  nforder *op = o->simplify();

  result->rtyp=nforder_type_id; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data
  Print("result is %lx\n", result->data);
  currRing->cf = (coeffs)result->data;

  return FALSE;
}





extern "C" int mod_init(SModulFunctions* psModulFunctions)
{
  nforder_Register();
  nforder_bb_setup();
  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),// the library name,
          "nfOrder",// the name for the singular interpreter
          FALSE,  // should not be static
          build_ring); // the C/C++ routine

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),// the library name,
          "pMaximalOrder",// the name for the singular interpreter
          FALSE,  // should not be static
          pMaximalOrder); // the C/C++ routine

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),// the library name,
          "oneStep",// the name for the singular interpreter
          FALSE,  // should not be static
          oneStep); // the C/C++ routine

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "calcdisc",
          FALSE, 
          _calcdisc); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "EltFromMat",
          FALSE, 
          elt_from_mat); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "NFOrderSimplify",
          FALSE, 
          nforder_simplify); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "EltNorm",
          FALSE, 
          elt_from_mat); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "EltTrace",
          FALSE, 
          elt_from_mat); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "EltRepMat",
          FALSE, 
          elt_from_mat); 


  module_help_main(
     (currPack->libname? currPack->libname: "NFOrder"),// the library name,
    "nforder: orders in number fields"); // the help string for the module
  return 1;
}
