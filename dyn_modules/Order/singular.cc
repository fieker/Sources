#include "kernel/mod2.h" // general settings/macros
#include"kernel/febase.h"  // for Print, WerrorS
#include"Singular/ipid.h" // for SModulFunctions, leftv
#include"libpolys/coeffs/numbers.h" // nRegister, coeffs.h
#include "libpolys/coeffs/coeffs.h"
#include"Singular/blackbox.h" // blackbox type
#include "nforder.h"
#include "libpolys/coeffs/bigintmat.h"

void dummy(coeffs) {}
static void WriteRing(const coeffs r, BOOLEAN details)
{
  printf("%s%d:%s: %s\n", __FILE__, __LINE__, __func__, "RING");
}

static char* CoeffString(const coeffs r)
{
  return "Ring";
}
static void WriteElt(number &a, const coeffs r)
{
  printf("%s%d:%s: %s\n", __FILE__, __LINE__, __func__, "RING-ELT");
}

static number EltCreate(long i, const coeffs r)
{
  return NULL;
}

static BOOLEAN order_cmp(coeffs, n_coeffType, void*)
{
  return FALSE;
}

static int nforder_type_id=0;
static n_coeffType nforder_type =n_CF;
BOOLEAN n_nfOrderInit(coeffs r,  void * parameter)
{
  puts("nfOrderInit called");
  printf("%s%d:%s: type is %d should be %d data %lx\n", __FILE__, __LINE__, __func__, getCoeffType(r), nforder_type, parameter);
  assume( getCoeffType(r) == nforder_type );
  r->nCoeffIsEqual=order_cmp;
  r->cfKillChar = NULL;
  r->cfCoeffString=CoeffString;
  r->cfCoeffWrite=WriteRing;
  r->cfWriteShort=WriteElt;
  r->cfInit = EltCreate;
  r->data = parameter;
  
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

  puts("build ring called");

  bigintmat **multtable = (bigintmat**)omAlloc(dimension*sizeof(bigintmat*));
  arg = arg->next;
  for (int i=0; i<dimension; i++) {
    multtable[i] = new bigintmat((bigintmat*)arg->Data());
    arg = arg->next;
  }
  nforder *o = new nforder(dimension, multtable, coeffs_BIGINT);
  result->rtyp=nforder_type_id; // set the result type
  result->data=(char*)nInitChar(nforder_type, o);// set the result data
  Print("result is %lx\n", result->data);

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
          "calcdisc",// the name for the singular interpreter
          FALSE,  // should not be static
          _calcdisc); // the C/C++ routine


  module_help_main(
     (currPack->libname? currPack->libname: "NFOrder"),// the library name,
    "nforder: orders in number fields"); // the help string for the module
  return 1;
}
