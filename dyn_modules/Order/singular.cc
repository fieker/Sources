#include "kernel/mod2.h" // general settings/macros
#include"kernel/febase.h"  // for Print, WerrorS
#include"Singular/ipid.h" // for SModulFunctions, leftv
#include"Singular/number2.h" // for SModulFunctions, leftv
#include"libpolys/coeffs/numbers.h" // nRegister, coeffs.h
#include "libpolys/coeffs/coeffs.h"
#include"Singular/blackbox.h" // blackbox type
#include "nforder.h"
#include "nforder_elt.h"
#include "nforder_ideal.h"
#include "libpolys/coeffs/bigintmat.h"

static int nforder_type_id=0;
n_coeffType nforder_type =n_unknown;

// coeffs stuff: -----------------------------------------------------------
static coeffs nforder_AE=NULL;
static void nforder_Register()
{
  puts("nforder_Register called");
  nforder_type=nRegister(n_unknown,n_nfOrderInit);
  nforder_AE=nInitChar(nforder_type,NULL);
}
// black box stuff: ---------------------------------------------------------
static void * nforder_ideal_Init(blackbox */*b*/)
{
  nforder_AE->ref++;
  return nforder_AE;
}
static char * nforder_ideal_String(blackbox *b, void *d)
{
  StringSetS("");
  if (d==NULL) StringAppendS("oo");
  else StringAppend("Ideal(%lx)",d);
  if (d) ((nforder_ideal *)d)->Write();
  return StringEndS();
}
static void * nforder_ideal_Copy(blackbox*b, void *d)
{  coeffs c=(coeffs)d; return new nforder_ideal((nforder_ideal*)d, 1);}

static BOOLEAN nforder_ideal_Assign(leftv l, leftv r)
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
static void nforder_ideal_destroy(blackbox *b, void *d)
{
  if (d!=NULL)
  {
    delete (nforder_ideal*)d;
  }
}

static BOOLEAN nforder_ideal_Op2(int op,leftv l, leftv r1, leftv r2)
{
  switch (op) {
    case '+':
      {
      nforder_ideal *I = (nforder_ideal*) r1->data,
                    *J = (nforder_ideal*) r2->data,
                    *H = nf_idAdd(I, J);
      l->rtyp = nforder_type_id;
      l->data = (void*)H;
      return FALSE;
      }
    default:
      return WrongOp("not implemented yet", op, r1);
  }
  return TRUE;
}
static BOOLEAN nforder_ideal_bb_setup()
{
  blackbox *b=(blackbox*)omAlloc0(sizeof(blackbox));
  // all undefined entries will be set to default in setBlackboxStuff
  // the default Print is quite useful,
  // all other are simply error messages
  b->blackbox_destroy=nforder_ideal_destroy;
  b->blackbox_String=nforder_ideal_String;
  //b->blackbox_Print=blackbox_default_Print;
  b->blackbox_Init=nforder_ideal_Init;
  b->blackbox_Copy=nforder_ideal_Copy;
  b->blackbox_Assign=nforder_ideal_Assign;
  //b->blackbox_Op1=blackbox_default_Op1;
  b->blackbox_Op2=nforder_ideal_Op2;
  //b->blackbox_Op3=blackbox_default_Op3;
  //b->blackbox_OpM=blackbox_default_OpM;
  nforder_type_id = setBlackboxStuff(b,"NFOrderIdeal");
  Print("setup: created a blackbox type [%d] '%s'",nforder_type_id, getBlackboxName(nforder_type_id));
  PrintLn();
  return FALSE; // ok, TRUE = error!
}

// module stuff: ------------------------------------------------------------

static BOOLEAN build_ring(leftv result, leftv arg)
{
  nforder *o;
  if (arg->Typ() == LIST_CMD) {
    lists L = (lists)arg->Data();
    int n = lSize(L)+1;
    bigintmat **multtable = (bigintmat**)omAlloc(n*sizeof(bigintmat*));
    for(int i=0; i<n; i++) {
      multtable[i] = (bigintmat*)(L->m[i].Data());
    }
    o = new nforder(n, multtable, nInitChar(n_Z, 0));
    omFree(multtable);
  } else {
    assume(arg->Typ() == INT_CMD);
    int dimension = (int)(long)arg->Data();

    bigintmat **multtable = (bigintmat**)omAlloc(dimension*sizeof(bigintmat*));
    arg = arg->next;
    for (int i=0; i<dimension; i++) {
      multtable[i] = new bigintmat((bigintmat*)arg->Data());
      arg = arg->next;
    }
    o = new nforder(dimension, multtable, nInitChar(n_Z, 0));
    for (int i=0; i<dimension; i++) {
      delete multtable[i];
    }
    omFree(multtable);
  }
  result->rtyp=CRING_CMD; // set the result type
  result->data=(char*)nInitChar(nforder_type, o);// set the result data

  return FALSE;
}

static BOOLEAN ideal_from_mat(leftv result, leftv arg)
{
  if (arg->Typ() != CRING_CMD) {
    WerrorS("1:usage: IdealFromMat(order, basis matrix)");
    return TRUE;
  }
  coeffs C = (coeffs)arg->Data();
  if (getCoeffType(C) != nforder_type) {
    WerrorS("2:usage: IdealFromMat(order, basis matrix)");
    return TRUE;
  }
  nforder *O = (nforder*) (C->data);
  arg = arg->next;
  if (arg->Typ()!=BIGINTMAT_CMD) {
    WerrorS("3:usage: IdealFromMat(order, basis matrix)");
    return TRUE;
  }
  bigintmat *b = (bigintmat*) arg->Data();
  result->rtyp = nforder_type_id;
  result->data = new nforder_ideal(b, C);
  return FALSE;
}


static BOOLEAN elt_from_mat(leftv result, leftv arg)
{
  number2 r = (number2)omAlloc(sizeof(struct snumber2));
  coeffs C = (coeffs)arg->Data();
  nforder *O = (nforder*) (C->data);
  arg = arg->next;
  bigintmat *b = (bigintmat*) arg->Data();
  result->rtyp = NUMBER2_CMD;
  r->n = (number)EltCreateMat(O, b);
  r->cf = nInitChar(nforder_type, O);
  result->data = r;
  return FALSE;
}

static BOOLEAN discriminant(leftv result, leftv arg)
{
  assume (arg->Typ()==CRING_CMD);
  coeffs c = (coeffs)arg->Data();
  assume (c->type == nforder_type);
  nforder * o = (nforder*)c->data;
  o->calcdisc();

  number2 tt = (number2)omAlloc(sizeof(struct snumber2));
  tt->n = o->getDisc();
  tt->cf = o->basecoeffs();
  result->rtyp = NUMBER2_CMD;
  result->data = tt;

  return FALSE;
}

static BOOLEAN pMaximalOrder(leftv result, leftv arg)
{
  assume (arg->Typ()==CRING_CMD);
  coeffs c = (coeffs)arg->Data();
  assume (c);
  assume (c->type == nforder_type);
  nforder * o = (nforder*)c->data;
  arg = arg->next;
  long p = (int)(long)arg->Data();
  number P = n_Init(p, o->basecoeffs());

  nforder *op = pmaximal(o, P);

  result->rtyp=CRING_CMD; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data
  assume(result->data);

  return FALSE;
}

static BOOLEAN oneStep(leftv result, leftv arg)
{
  assume (arg->Typ()==CRING_CMD);
  coeffs c = (coeffs)arg->Data();
  assume (c->type == nforder_type);
  nforder * o = (nforder*)c->data;
  arg = arg->next;
  long p = (int)(long)arg->Data();
  number P = n_Init(p, o->basecoeffs());

  nforder *op = onestep(o, P, o->basecoeffs());

  result->rtyp=CRING_CMD; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data

  return FALSE;
}

static BOOLEAN nforder_simplify(leftv result, leftv arg)
{
  assume (arg->Typ()==CRING_CMD);
  coeffs c = (coeffs)arg->Data();
  assume (c->type = nforder_type);
  nforder * o = (nforder*)c->data;

  nforder *op = o->simplify();

  result->rtyp=CRING_CMD; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data

  return FALSE;
}

static BOOLEAN eltTrace(leftv result, leftv arg)
{
  assume (arg->Typ()==NUMBER2_CMD);
  number2 a = (number2) arg->Data();
  coeffs  c = a->cf;
  bigintmat * aa = (bigintmat*)a->n;
  assume (c->type == nforder_type);
  nforder * o = (nforder*)c->data;
  number t = o->elTrace(aa);
  number2 tt = (number2)omAlloc(sizeof(struct snumber2));
  tt->n = t;
  tt->cf = o->basecoeffs();
  result->rtyp = NUMBER2_CMD;
  result->data = tt;
  return FALSE;
}

static BOOLEAN eltNorm(leftv result, leftv arg)
{
  assume (arg->Typ()==NUMBER2_CMD);
  number2 a = (number2) arg->Data();
  coeffs  c = a->cf;
  bigintmat * aa = (bigintmat*)a->n;
  assume (c->type == nforder_type);
  nforder * o = (nforder*)c->data;
  number t = o->elNorm(aa);
  number2 tt = (number2)omAlloc(sizeof(struct snumber2));
  tt->n = t;
  tt->cf = o->basecoeffs();
  result->rtyp = NUMBER2_CMD;
  result->data = tt;
  return FALSE;
}

static BOOLEAN eltRepMat(leftv result, leftv arg)
{
  assume (arg->Typ()==NUMBER2_CMD);
  number2 a = (number2) arg->Data();
  coeffs  c = a->cf;
  bigintmat * aa = (bigintmat*)a->n;
  assume (c->type == nforder_type);
  nforder * o = (nforder*)c->data;
  bigintmat* t = o->elRepMat(aa);
  result->rtyp = BIGINTMAT_CMD;
  result->data = t;
  return FALSE;
}

static BOOLEAN smithtest(leftv result, leftv arg)
{
  assume (arg->Typ()==BIGINTMAT_CMD);
  bigintmat *a = (bigintmat *) arg->Data();
  arg = arg->next;

  long p = (int)(long)arg->Data();
  number P = n_Init(p, a->basecoeffs());

  bigintmat * A, *B;
  diagonalForm(a, &A, &B);
 

  result->rtyp = NONE;
  return FALSE;
}


extern "C" int mod_init(SModulFunctions* psModulFunctions)
{
  nforder_Register();
  nforder_ideal_bb_setup();
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
          "Discriminant",
          FALSE, 
          discriminant); 

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
          eltNorm); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "EltTrace",
          FALSE, 
          eltTrace); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "EltRepMat",
          FALSE, 
          eltRepMat); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "SmithTest",
          FALSE, 
          smithtest); 

  psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          "IdealFromMat",
          FALSE, 
          ideal_from_mat); 

  module_help_main(
     (currPack->libname? currPack->libname: "NFOrder"),// the library name,
    "nforder: orders in number fields"); // the help string for the module
  return 1;
}
