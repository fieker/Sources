#include "kernel/mod2.h" // general settings/macros
#include "Singular/mod_lib.h"
#include"Singular/ipid.h" // for SModulFunctions, leftv
#include"Singular/number2.h" // for SModulFunctions, leftv
#include"libpolys/coeffs/numbers.h" // nRegister, coeffs.h
#include "libpolys/coeffs/coeffs.h"
#include"Singular/blackbox.h" // blackbox type
#include <Singular/libsingular.h> // for iiLibCmd

#include "nforder.h"
#include "nforder_elt.h"
#include "nforder_ideal.h"
#include "lattice.h"
#include "libpolys/coeffs/bigintmat.h"

#ifdef SINGULAR_4_1
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
  if (d) ((nforder_ideal *)d)->Write();
  else StringAppendS("o not defined o");
  return StringEndS();
}
static void * nforder_ideal_Copy(blackbox* /*b*/, void *d)
{ return new nforder_ideal((nforder_ideal*)d, 1);}

static BOOLEAN nforder_ideal_Assign(leftv l, leftv r)
{
  if (l->Typ()==r->Typ())
  {
    if (l->rtyp==IDHDL)
    {
      IDDATA((idhdl)l->data)=(char *)nforder_ideal_Copy((blackbox*)NULL, r->data);
    }
    else
    {
      l->data=(char *)nforder_ideal_Copy((blackbox*)NULL, r->data);
    }
    return FALSE;
  }
  return TRUE;
}
static void nforder_ideal_destroy(blackbox * /*b*/, void *d)
{
  if (d!=NULL)
  {
    delete (nforder_ideal*)d;
  }
}

BOOLEAN checkArgumentIsOrder(leftv arg, nforder * cmp, nforder ** result)
{
  if (arg->Typ() != CRING_CMD) return FALSE;
  coeffs R = (coeffs) arg->Data();
  if (getCoeffType(R) != nforder_type) return FALSE;
  nforder * O = (nforder*) R->data;
  if (cmp && cmp != O) return FALSE;
  *result = O;
  return TRUE;
}

BOOLEAN checkArgumentIsBigintmat(leftv arg, coeffs r, bigintmat ** result)
{
  if (arg->Typ() != BIGINTMAT_CMD) return FALSE;
  bigintmat * b = (bigintmat*) arg->Data();
  if (r && b->basecoeffs() != r) return FALSE;
  *result = b;
  return TRUE;
}

BOOLEAN checkArgumentIsNumber2(leftv arg, coeffs r, number2 * result)
{
  if (arg->Typ() != CNUMBER_CMD) return FALSE;
  number2 b = (number2) arg->Data();
  if (r && b->cf != r) return FALSE;
  *result = b;
  return TRUE;
}


BOOLEAN checkArgumentIsNFOrderIdeal(leftv arg, coeffs r, nforder_ideal ** result)
{
  if (arg->Typ() != nforder_type_id) return FALSE;
  *result = (nforder_ideal *) arg->Data();
  if (r && (*result)->order() != r) return FALSE;
  return TRUE;
}

BOOLEAN checkArgumentIsInt(leftv arg, int* result)
{
  if (arg->Typ() != INT_CMD) return FALSE;
  *result = (long) arg->Data();
  return TRUE;
}

BOOLEAN checkArgumentIsBigint(leftv arg, number* result)
{
  switch (arg->Typ()) {
    case BIGINT_CMD:
      *result = (number)arg->Data();
      return TRUE;
      break;
    case NUMBER_CMD:
      if (currRing->cf == coeffs_BIGINT &&
          getCoeffType(coeffs_BIGINT) == n_Z) {
        *result = (number)arg->Data();
        return TRUE;
      } else
        return FALSE;
      break;
    case CNUMBER_CMD: 
      {
        number2 n = (number2)arg->Data();
        if (getCoeffType(n->cf) == n_Z) {
          *result = n->n;
          return TRUE;
        } 
        return FALSE;
        break;
      }
    default:
      return FALSE;
  }
}

static BOOLEAN nforder_ideal_Op2(int op,leftv l, leftv r1, leftv r2)
{
  Print("Types are %d %d\n", r1->Typ(), r2->Typ());
  number2 e;
  int f;
  nforder_ideal *I, *J, *H;
  switch (op) {
    case '+':
      {
      if (!checkArgumentIsNFOrderIdeal(r1, NULL, &I))
        return TRUE;
      if (!checkArgumentIsNFOrderIdeal(r2, I->order(), &J))
        return TRUE;
      H = nf_idAdd(I, J);
      break;
      }
    case '*':
      {
      if (!checkArgumentIsNFOrderIdeal(r1, NULL, &I)) {
        leftv r = r1;
        r1 = r2;
        r2 = r; //at least ONE argument has to be an ideal
      }
      if (!checkArgumentIsNFOrderIdeal(r1, NULL, &I))
        return TRUE;
      if (checkArgumentIsNFOrderIdeal(r2, I->order(), &J)) {
        H = nf_idMult(I, J);
      } else if (checkArgumentIsNumber2(r2, I->order(), &e)) {
        H = nf_idMult(I, e->n);
      } else if (checkArgumentIsInt(r2, &f)) {
        H = nf_idMult(I, f);
      } else
        return TRUE;
      break;
      }
    case '^':
      {
        if (!checkArgumentIsNFOrderIdeal(r1, NULL, &I))
          return TRUE;
        if (!checkArgumentIsInt(r2, &f))
          return TRUE;
        H = nf_idPower(I, f);
        break;
      }
    default:
      return TRUE;
  }
  l->rtyp = nforder_type_id;
  l->data = (void*)H;
  return FALSE;
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

// lattice type: ------------------------------------------------------------

static int lattice_id=1; //NOTE 1 works?

static void * lattice_Init(blackbox */*b*/)
{
  void * l = omAlloc(sizeof(lattice));
//   Print("create at %lx\n",(unsigned long)l);
  return (void*) l;
}

static void lattice_destroy(blackbox * /*b*/, void *d)
{
  if (d!=NULL)
  {
//     Print("destroy %x at %lx\n",*((int*)d),(unsigned long)d);
    delete (lattice*)d;
  }
}

static void * lattice_Copy(blackbox* /*b*/, void *d)
{ 
  lattice * ll = (lattice*) omAlloc(sizeof(lattice));
  *ll = *((lattice*)d);
  return ll;   //NOTE
}


static char * lattice_String(blackbox */*b*/, void *d)
{
  StringSetS("");
  if(d){
    ((lattice *)d)->Write();
  } else {
    StringAppendS("o not defined o");
  }
  return StringEndS();
}
static BOOLEAN lattice_Assign(leftv l, leftv r)
{
  if (l->Typ()==r->Typ()) // assignment of same type
  {
    if (l->rtyp==IDHDL) {// assign to variable    
      omFree(IDDATA((idhdl)l->data));
      IDDATA((idhdl)l->data)=(char *)lattice_Copy((blackbox*)NULL, r->data);
    } else { // assign to a part of a structure
      leftv ll = l->LData();
      if(ll == NULL) {
          return TRUE; // out of array bounds or similiar
      }
      omFree(ll->data);
      l->data=(char*)lattice_Copy((blackbox*)NULL, r->data);
    }
//   } else if (r->Typ() == BIGINTMAT_CMD) { // use matrix as lattice base
//     bigintmat * basis = (bigintmat*)r->Data();
//     lattice * lat = new lattice(basis,false);
// //     omFree(IDDATA((idhdl)l->data));
//     l->data=(char*) lat;
  } else {
    Werror("assign Type(%d) = Type(%d) not implemented",l->Typ(),r->Typ());
    return TRUE;
  }
  return FALSE;
}

static BOOLEAN lattice_bb_setup()
{
  blackbox * b = (blackbox*) omAlloc0(sizeof(blackbox));
  b->blackbox_Init = lattice_Init;
  b->blackbox_destroy = lattice_destroy;
  b->blackbox_Copy = lattice_Copy;
  b->blackbox_String = lattice_String;
//   b->blackbox_Print = lattice_Print;
  b->blackbox_Assign = lattice_Assign;
  lattice_id = setBlackboxStuff(b,"lattice");    

  Print("setup: created a blackbox type [%d] '%s'",lattice_id, getBlackboxName(lattice_id));
  PrintLn();
  return FALSE;
}



// module stuff: ------------------------------------------------------------

BOOLEAN checkBigintmatDim(bigintmat * b, int r, int c)
{
  if (b->rows() != r) return FALSE;
  if (b->cols() != c) return FALSE;
  return TRUE;
}

#define returnNumber(_res, _n, _R) \
  do {                                                          \
    number2 _r = (number2)omAlloc(sizeof(struct snumber2));     \
    _r->n = _n;                                                 \
    _r->cf = _R;                                                \
    _res->rtyp = CNUMBER_CMD;                                   \
    _res->data = _r;                                            \
  } while (0)
    

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
  nforder * O;
  if (!checkArgumentIsOrder(arg, NULL, &O)) {
    WerrorS("usage: IdealFromMat(order, basis matrix)");
    return TRUE;
  }
  arg = arg->next;
  bigintmat *b;
  if (!checkArgumentIsBigintmat(arg, O->basecoeffs(), &b)) {
    WerrorS("3:usage: IdealFromMat(order, basis matrix)");
    return TRUE;
  }
  result->rtyp = nforder_type_id;
  result->data = new nforder_ideal(b, nInitChar(nforder_type, O));
  return FALSE;
}


static BOOLEAN elt_from_mat(leftv result, leftv arg)
{
  nforder * O;
  if (!checkArgumentIsOrder(arg, NULL, &O)) {
    WerrorS("usage: EltFromMat(order, matrix)");
    return TRUE;
  }
  arg = arg->next;
  bigintmat *b;
  if (!checkArgumentIsBigintmat(arg, O->basecoeffs(), &b)) {
    WerrorS("2:usage: EltFromMat(order, matrix)");
    return TRUE;
  }
  returnNumber(result, (number)EltCreateMat(O, b), nInitChar(nforder_type, O));
  return FALSE;
}

static BOOLEAN discriminant(leftv result, leftv arg)
{
  nforder * O;
  if (!checkArgumentIsOrder(arg, NULL, &O)) {
    WerrorS("usage: Discriminant(order)");
    return TRUE;
  }
  O->calcdisc();

  returnNumber(result, O->getDisc(), O->basecoeffs());
  return FALSE;
}

static BOOLEAN pMaximalOrder(leftv result, leftv arg)
{
  nforder * o;
  if (!checkArgumentIsOrder(arg, NULL, &o)) {
    WerrorS("usage: pMaximalOrder(order, int)");
    return TRUE;
  }
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
  nforder * o;
  if (!checkArgumentIsOrder(arg, NULL, &o)) {
    WerrorS("usage: NFOrderSimplify(order)");
    return TRUE;
  }
  nforder *op = o->simplify();

  result->rtyp=CRING_CMD; // set the result type
  result->data=(char*)nInitChar(nforder_type, op);// set the result data

  return FALSE;
}

static BOOLEAN eltTrace(leftv result, leftv arg)
{
  number2 a;
  if (!checkArgumentIsNumber2(arg, NULL, &a)) {
    WerrorS("EltTrace(elt)");
    return TRUE;
  }
  coeffs  c = a->cf;
  if (getCoeffType(c) != nforder_type) {
    WerrorS("EltTrace(elt in order)");
    return TRUE;
  }
  bigintmat * aa = (bigintmat*)a->n;
  nforder * o = (nforder*)c->data;
  number t = o->elTrace(aa);
  returnNumber(result, t, o->basecoeffs());
  return FALSE;
}

static BOOLEAN eltNorm(leftv result, leftv arg)
{
  number2 a;
  if (!checkArgumentIsNumber2(arg, NULL, &a)) {
    WerrorS("EltNorm(elt)");
    return TRUE;
  }
  coeffs  c = a->cf;
  if (getCoeffType(c) != nforder_type) {
    WerrorS("EltNorm(elt in order)");
    return TRUE;
  }
  bigintmat * aa = (bigintmat*)a->n;
  nforder * o = (nforder*)c->data;
  number t = o->elNorm(aa);
  returnNumber(result, t, o->basecoeffs());
  return FALSE;
}

static BOOLEAN eltRepMat(leftv result, leftv arg)
{
  assume (arg->Typ()==CNUMBER_CMD);
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

static BOOLEAN bimToCurrRing(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != BIGINTMAT_CMD)) 
  {
    WerrorS("usage: bimToCurrRing(bigintmat)");
    return TRUE;
  }

  bigintmat * in = (bigintmat *) arg->Data();  
  bigintmat * out = bimChangeCoeff(in,currRing->cf);

  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) out;
  return FALSE;
}

static BOOLEAN latticeFromBasis(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != BIGINTMAT_CMD)) 
  {
    WerrorS("usage: latticeFromBasis(bigintmat)");
    return TRUE;
  }

  bigintmat * a = (bigintmat *) arg->Data();  
  lattice * l = new lattice(a,false);

  result->rtyp = lattice_id;
  result->data = (void*) l;
  return FALSE;
}
static BOOLEAN latticeFromGramMatrix(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != BIGINTMAT_CMD)) 
  {
    WerrorS("usage: latticeFromGramMatrix(bigintmat)");
    return TRUE;
  }

  bigintmat * a = (bigintmat *) arg->Data();  
  lattice * l = new lattice(a,true);

  result->rtyp = lattice_id;
  result->data = (void*) l;
  return FALSE;
}

static BOOLEAN LLL(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: LLL(lattice,[number])");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();
  
  number c;
  coeffs coef = currRing->cf;
  
  arg = arg->next;
  
  if(arg == NULL) 
  {
    number three = n_Init(3, coef);
    number four  = n_Init(4, coef);
    c = n_Div(three,four,coef);
    n_Delete(&three,coef);
    n_Delete(&four,coef);
  } 
  else if(arg->Typ() != NUMBER_CMD) 
  {
    WerrorS("usage: LLL(lattice,[number])");
    return TRUE;
  } else {
    c = n_Copy((number) arg->Data(),coef);
  }
  
  l->LLL(c,coef,true,false,true);
  n_Delete(&c,coef);
  
  result->rtyp = NONE;
  return FALSE;
}

static BOOLEAN getLatticeElement(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: getLatticeElement(lattice, bigintmat)");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();
  arg = arg->next;
  if( (arg == NULL) 
    ||(arg->Typ() != BIGINTMAT_CMD)) 
  {
    WerrorS("usage: enumerateAll(lattice, number)");
    return TRUE;
  }
  bigintmat * in = ((bigintmat *)arg->Data());
  bigintmat * enumeration = l->get_lattice_element(in);
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) enumeration;
  delete in;
  return FALSE;
}
static BOOLEAN getBasis(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: getBasis(lattice)");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();

  bigintmat * basis = l->get_basis();
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) basis;
  return FALSE;
}

static BOOLEAN getReducedBasis(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: getReducedBasis(lattice)");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();

  bigintmat * reduced = l->get_reduced_basis();
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) reduced;
  return FALSE;
}

static BOOLEAN getTransformationMatrix(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: getTransformationMatrix(lattice)");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();

  bigintmat * H = l->get_transformation_matrix();
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) H;
  return FALSE;
}
static BOOLEAN getGramMatrix(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: getGramMatrix(lattice)");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();

  bigintmat * gram = l->get_gram_matrix();
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) gram;
  return FALSE;
}

static BOOLEAN enumerateAll(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: enumerateAll(lattice, number)");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();
  arg = arg->next;
  if( (arg == NULL) 
    ||(arg->Typ() != NUMBER_CMD)) 
  {
    WerrorS("usage: enumerateAll(lattice, number)");
    return TRUE;
  }
  number c = ((number)arg->Data());
  bigintmat * enumeration = l->enumerate_all(c);
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) enumeration;
  return FALSE;
}

static BOOLEAN enumerateNext(leftv result, leftv arg)
{ 
  if( (arg == NULL) 
    ||(arg->Typ() != lattice_id)) 
  {
    WerrorS("usage: enumerateNext(lattice [, number] [,bigintmat])");
    return TRUE;
  }
  lattice * l = (lattice*) arg->Data();
  arg = arg->next;
  bigintmat * enumeration = NULL;
  
  
  if( (arg == NULL) || ((arg->Typ() != NUMBER_CMD) && (arg->Typ() != BIGINTMAT_CMD)) ) 
  {
    enumeration = l->enumerate_next();
  } else {
    if(arg->Typ() == NUMBER_CMD)
    {
      number c = ((number)arg->Data());
      arg = arg->next;
      if( (arg == NULL) )
      {
        enumeration = l->enumerate_next(c);
      } else {
        if(arg->Typ() == BIGINTMAT_CMD)
        {
          bigintmat * in = (bigintmat *) arg->Data();
          enumeration = l->enumerate_next(c,in);
        } else {
          WerrorS("usage: enumerateNext(lattice [, number] [,bigintmat])");
          return TRUE;
        }
      }
    } else {
      if(arg->Typ() == BIGINTMAT_CMD)
      {
        bigintmat * in = (bigintmat *) arg->Data();
        enumeration = l->enumerate_next(in);
      } else {
        WerrorS("usage: enumerateNext(lattice [, number] [,bigintmat])");
        return TRUE;
      }
    }
  }
  if(enumeration == NULL)
  {
    bigintmat * basis =l->get_basis();
    enumeration = new bigintmat(basis->cols(),1,basis->basecoeffs());
    delete basis;
  }
  
  result->rtyp = BIGINTMAT_CMD;
  result->data = (void*) enumeration;
  return FALSE;
}

static BOOLEAN get_same_field_poly(leftv result, leftv arg)
{
  if( (arg == NULL) 
    ||(arg->Typ() != POLY_CMD)) 
  {
    WerrorS("usage: sameFieldPoly( poly )");
    return TRUE;
  }
  poly in =((poly)arg->Data());
  poly out = get_nice_poly(in);
  result->rtyp = POLY_CMD;
  result->data = (void*) out;
  return FALSE;
}
static BOOLEAN t2_norm(leftv result, leftv arg)
{
  if( (arg == NULL) 
    ||(arg->Typ() != POLY_CMD)) 
  {
    WerrorS("usage: t2norm(poly, int)");
    return TRUE;
  }
  poly in =((poly)arg->Data());
  arg = arg->next;
  if( (arg == NULL) 
    ||(arg->Typ() != INT_CMD)) 
  {
    WerrorS("usage: t2norm(poly, int)");
    return TRUE;
  }
  int prec = (int)(long) arg->Data();
  number out = t2norm(in,currRing,currRing->cf,prec);
  result->rtyp = NUMBER_CMD;
  result->data = (void*) out;
  return FALSE;
}

typedef struct {
  const char * name;
  BOOLEAN (*func)(leftv, leftv);
} cmd;

static cmd cmd_arr[] = {
{ "nfOrder", build_ring},
{ "pMaximalOrder", pMaximalOrder},
{ "OneStep", oneStep},
{ "Discriminant", discriminant}, 
{ "EltFromMat", elt_from_mat}, 
{ "NFOrderSimplify", nforder_simplify}, 
{ "EltNorm", eltNorm}, 
{ "EltTrace", eltTrace}, 
{ "EltRepMat", eltRepMat}, 
{ "SmithTest", smithtest}, 
{ "IdealFromMat", ideal_from_mat}, 
{ "bimToCurrRing", bimToCurrRing},
{ "latticeFromBasis", latticeFromBasis},
{ "latticeFromGramMatrix", latticeFromGramMatrix},
{ "LLL", LLL},
{ "getLatticeElement", getLatticeElement},
{ "getBasis", getBasis},
{ "getReducedBasis", getReducedBasis},
{ "getTransformationMatrix", getTransformationMatrix},
{ "getGramMatrix", getGramMatrix}, 
{ "enumerateAll", enumerateAll},
{ "enumerateNext", enumerateNext},
{ "sameFieldPoly", get_same_field_poly},
{ "t2norm", t2_norm},
{ NULL, NULL}
};

extern "C" int SI_MOD_INIT(Order)(SModulFunctions* psModulFunctions)
{
  //load nforder.lib for additional procedures
  //NOTE: is this the correct way to do it?
  iiLibCmd(omStrDup("nforder.lib"), TRUE,TRUE,TRUE);

  nforder_Register();

  nforder_ideal_bb_setup();
  lattice_bb_setup();


  for(int i=0; cmd_arr[i].name; i++) {
    psModulFunctions->iiAddCproc(
          (currPack->libname? currPack->libname: ""),
          cmd_arr[i].name,
          FALSE, 
          cmd_arr[i].func);
  }

  module_help_main(
     (currPack->libname? currPack->libname: "NFOrder"),// the library name,
    "nforder: orders in number fields"); // the help string for the module
  return MAX_TOK;
}
#endif
