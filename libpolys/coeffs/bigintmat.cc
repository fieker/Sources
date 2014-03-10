/*****************************************
 *  Computer Algebra System SINGULAR      *
 *****************************************/
/*
 *  * ABSTRACT: class bigintmat: matrices of numbers.
 *   * a few functinos might be limited to bigint or euclidean rings.
 *    */

#ifdef HAVE_CONFIG_H
#include "libpolysconfig.h"
#endif /* HAVE_CONFIG_H */


#include "bigintmat.h"
#include <misc/intvec.h>

#include <math.h>
#include <string.h>


//#define BIMATELEM(M,I,J) (M)[ (M).index(I,J) ]

bigintmat * bigintmat::transpose()
{
  bigintmat * t = new bigintmat(col, row, basecoeffs());
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      t->set(j, i, BIMATELEM(*this,i,j));
    }
  }
  return t;
}

// Beginnt bei [0]
void bigintmat::set(int i, number n, const coeffs C)
{
  assume (C == NULL || C == basecoeffs());

  rawset(i, n_Copy(n, basecoeffs()), basecoeffs());
}

// Beginnt bei [1,1]
void bigintmat::set(int i, int j, number n, const coeffs C)
{
  assume (C == NULL || C == basecoeffs());
  assume (i >= 0 && j >= 0);
  assume (i <= rows() && j <= cols());
  set(index(i, j), n, C);
}

number bigintmat::get(int i) const
{
  assume (i >= 0);
  assume (i<rows()*cols());

  return n_Copy(v[i], basecoeffs());
}

number bigintmat::view(int i) const
{
  assume (i >= 0);
  assume (i<rows()*cols());

  return v[i];
}

number bigintmat::get(int i, int j) const
{
  assume (i >= 0 && j >= 0);
  assume (i <= rows() && j <= cols());

  return get(index(i, j));
}

number bigintmat::view(int i, int j) const
{
  assume (i >= 0 && j >= 0);
  assume (i <= rows() && j <= cols());

  return view(index(i, j));
}
// Ueberladener *=-Operator (für int und bigint)
// Frage hier: *= verwenden oder lieber = und * einzeln?
void bigintmat::operator*=(int intop)
{
  number iop = n_Init(intop, basecoeffs());

  inpMult(iop, basecoeffs());

  n_Delete(&iop, basecoeffs());
}

void bigintmat::inpMult(number bintop, const coeffs C)
{
  assume (C == NULL || C == basecoeffs());

  const int l = rows() * cols();

  for (int i=0; i < l; i++)
    n_InpMult(v[i], bintop, basecoeffs());
}

// Stimmen Parameter?
// Welche der beiden Methoden?
// Oder lieber eine comp-Funktion?

bool operator==(const bigintmat & lhr, const bigintmat & rhr)
{
  if (&lhr == &rhr) { return true; }
  if (lhr.cols() != rhr.cols()) { return false; }
  if (lhr.rows() != rhr.rows()) { return false; }
  if (lhr.basecoeffs() != rhr.basecoeffs()) { return false; }

  const int l = (lhr.rows())*(lhr.cols());

  for (int i=0; i < l; i++)
  {
    if (!n_Equal(lhr[i], rhr[i], lhr.basecoeffs())) { return false; }
  }

  return true;
}

bool operator!=(const bigintmat & lhr, const bigintmat & rhr)
{
  return !(lhr==rhr);
}

// Matrix-Add/-Sub/-Mult so oder mit operator+/-/* ?
bigintmat * bimAdd(bigintmat * a, bigintmat * b)
{
  if (a->cols() != b->cols()) return NULL;
  if (a->rows() != b->rows()) return NULL;
  if (a->basecoeffs() != b->basecoeffs()) { return NULL; }

  const coeffs basecoeffs = a->basecoeffs();

  int i;

  bigintmat * bim = new bigintmat(a->rows(), a->cols(), basecoeffs);

  for (i=a->rows()*a->cols()-1;i>=0; i--)
    bim->rawset(i, n_Add((*a)[i], (*b)[i], basecoeffs), basecoeffs);

  return bim;
}
bigintmat * bimAdd(bigintmat * a, int b)
{

  const int mn = a->rows()*a->cols();

  const coeffs basecoeffs = a->basecoeffs();
  number bb=n_Init(b,basecoeffs);

  int i;

  bigintmat * bim = new bigintmat(a->rows(),a->cols() , basecoeffs);

  for (i=0; i<mn; i++)
    bim->rawset(i, n_Add((*a)[i], bb, basecoeffs), basecoeffs);

  n_Delete(&bb,basecoeffs);
  return bim;
}

bigintmat * bimSub(bigintmat * a, bigintmat * b)
{
  if (a->cols() != b->cols()) return NULL;
  if (a->rows() != b->rows()) return NULL;
  if (a->basecoeffs() != b->basecoeffs()) { return NULL; }

  const coeffs basecoeffs = a->basecoeffs();

  int i;

  bigintmat * bim = new bigintmat(a->rows(), a->cols(), basecoeffs);

  for (i=a->rows()*a->cols()-1;i>=0; i--)
    bim->rawset(i, n_Sub((*a)[i], (*b)[i], basecoeffs), basecoeffs);

  return bim;
}

bigintmat * bimSub(bigintmat * a, int b)
{

  const int mn = a->rows()*a->cols();

  const coeffs basecoeffs = a->basecoeffs();
  number bb=n_Init(b,basecoeffs);

  int i;

  bigintmat * bim = new bigintmat(a->rows(),a->cols() , basecoeffs);

  for (i=0; i<mn; i++)
    bim->rawset(i, n_Sub((*a)[i], bb, basecoeffs), basecoeffs);

  n_Delete(&bb,basecoeffs);
  return bim;
}

bigintmat * bimMult(bigintmat * a, bigintmat * b)
{
  const int ca = a->cols();
  const int cb = b->cols();

  const int ra = a->rows();
  const int rb = b->rows();

  if (ca != rb)
  {
#ifndef SING_NDEBUG
    Werror("wrong bigintmat sizes at multiplication a * b: acols: %d != brows: %d\n", ca, rb);
#endif
    return NULL;
  }

  assume (ca == rb);

  if (a->basecoeffs() != b->basecoeffs()) { return NULL; }

  const coeffs basecoeffs = a->basecoeffs();

  int i, j, k;

  number sum;

  bigintmat * bim = new bigintmat(ra, cb, basecoeffs);

  for (i=1; i<=ra; i++)
    for (j=1; j<=cb; j++)
    {
      sum = n_Init(0, basecoeffs);

      for (k=1; k<=ca; k++)
      {
        number prod = n_Mult( BIMATELEM(*a, i, k), BIMATELEM(*b, k, j), basecoeffs);

        number sum2 = n_Add(sum, prod, basecoeffs); // no inplace add :(

        n_Delete(&sum, basecoeffs); n_Delete(&prod, basecoeffs);

        sum = sum2;
      }
      bim->rawset(i, j, sum, basecoeffs);
    }
  return bim;
}

bigintmat * bimMult(bigintmat * a, int b)
{

  const int mn = a->rows()*a->cols();

  const coeffs basecoeffs = a->basecoeffs();
  number bb=n_Init(b,basecoeffs);

  int i;

  bigintmat * bim = new bigintmat(a->rows(),a->cols() , basecoeffs);

  for (i=0; i<mn; i++)
    bim->rawset(i, n_Mult((*a)[i], bb, basecoeffs), basecoeffs);

  n_Delete(&bb,basecoeffs);
  return bim;
}

bigintmat * bimMult(bigintmat * a, number b, const coeffs cf)
{
  if (cf!=a->basecoeffs()) return NULL;

  const int mn = a->rows()*a->cols();

  const coeffs basecoeffs = a->basecoeffs();

  int i;

  bigintmat * bim = new bigintmat(a->rows(),a->cols() , basecoeffs);

  for (i=0; i<mn; i++)
    bim->rawset(i, n_Mult((*a)[i], b, basecoeffs), basecoeffs);

  return bim;
}

// ----------------------------------------------------------------- //
// Korrekt?

intvec * bim2iv(bigintmat * b)
{
  intvec * iv = new intvec(b->rows(), b->cols(), 0);
  for (int i=0; i<(b->rows())*(b->cols()); i++)
    (*iv)[i] = n_Int((*b)[i], b->basecoeffs()); // Geht das so?
  return iv;
}

bigintmat * iv2bim(intvec * b, const coeffs C)
{
  const int l = (b->rows())*(b->cols());
  bigintmat * bim = new bigintmat(b->rows(), b->cols(), C);

  for (int i=0; i < l; i++)
    bim->rawset(i, n_Init((*b)[i], C), C);

  return bim;
}

// ----------------------------------------------------------------- //

int bigintmat::compare(const bigintmat* op) const
{
  assume (basecoeffs() == op->basecoeffs() );

#ifndef SING_NDEBUG
  if (basecoeffs() != op->basecoeffs() )
    WerrorS("wrong bigintmat comparison: different basecoeffs!\n");
#endif

  if ((col!=1) ||(op->cols()!=1))
  {
    if((col!=op->cols())
       || (row!=op->rows()))
      return -2;
  }

  int i;
  for (i=0; i<si_min(row*col,op->rows()*op->cols()); i++)
  {
    if ( n_Greater(v[i], (*op)[i], basecoeffs()) )
      return 1;
    else if (! n_Equal(v[i], (*op)[i], basecoeffs()))
      return -1;
  }

  for (; i<row; i++)
  {
    if ( n_GreaterZero(v[i], basecoeffs()) )
      return 1;
    else if (! n_IsZero(v[i], basecoeffs()) )
      return -1;
  }
  for (; i<op->rows(); i++)
  {
    if ( n_GreaterZero((*op)[i], basecoeffs()) )
      return -1;
    else if (! n_IsZero((*op)[i], basecoeffs()) )
      return 1;
  }
  return 0;
}


bigintmat * bimCopy(const bigintmat * b)
{
  if (b == NULL)
    return NULL;

  return new bigintmat(b);
}

void bigintmat::Write()
{
  int n = cols(), m=rows();

  StringAppendS("[ ");
  for(int i=1; i<= m; i++) {
    StringAppendS("[ ");
    for(int j=1; j< n; j++) {
      n_Write(v[(i-1)*n+j-1], basecoeffs());
      StringAppendS(", ");
    }
    n_Write(v[i*n-1], basecoeffs());
    StringAppendS(" ]");
    if (i<m) {
      StringAppendS(", ");
    }
  }
  StringAppendS(" ] ");
}

char* bigintmat::String()
{
  StringSetS("");
  Write();
  return StringEndS();
}

void bigintmat::Print()
{
  StringSetS("");
  Write();
  ::Print("%s\n", StringEndS());
}


char* bigintmat::StringAsPrinted()
{
  if ((col==0) || (row==0)) 
    return NULL;
  else
  {
    int * colwid = getwid(80);
    if (colwid == NULL)
    {
      WerrorS("not enough space to print bigintmat");
      WerrorS("try string(...) for a unformatted output");
      return NULL;
    }
    char * ps;
    int slength = 0;
    for (int j=0; j<col; j++)
      slength += colwid[j]*row;
    slength += col*row+row;
    ps = (char*) omAlloc0(sizeof(char)*(slength));
    int pos = 0;
    for (int i=0; i<col*row; i++)
    {
      StringSetS("");
      n_Write(v[i], basecoeffs());
      char * ts = StringEndS();
      const int _nl = strlen(ts);
      int cj = i%col;
      if (_nl > colwid[cj])
      {
        StringSetS("");
        int ci = i/col;
        StringAppend("[%d,%d]", ci+1, cj+1);
        char * ph = StringEndS();
        int phl = strlen(ph);
        if (phl > colwid[cj])
        {
          for (int j=0; j<colwid[cj]-1; j++)
            ps[pos+j] = ' ';
          ps[pos+colwid[cj]-1] = '*';
        }
        else
        {
          for (int j=0; j<colwid[cj]-phl; j++)
            ps[pos+j] = ' ';
          for (int j=0; j<phl; j++)
            ps[pos+colwid[cj]-phl+j] = ph[j];
        }
        omFree(ph);
      }
      else  // Mit Leerzeichen auffüllen und zahl reinschreiben
      {
        for (int j=0; j<colwid[cj]-_nl; j++)
          ps[pos+j] = ' ';
        for (int j=0; j<_nl; j++)
          ps[pos+colwid[cj]-_nl+j] = ts[j];
      }
      // ", " und (evtl) "\n" einfügen
      if ((i+1)%col == 0)
      {
        if (i != col*row-1)
        {
          ps[pos+colwid[cj]] = ',';
          ps[pos+colwid[cj]+1] = '\n';
          pos += colwid[cj]+2;
        }
      }
      else
      {
        ps[pos+colwid[cj]] = ',';
        pos += colwid[cj]+1;
      }

      omFree(ts);  // Hier ts zerstören
    }
    return(ps);
    // omFree(ps);
  }
}

int intArrSum(int * a, int length)
{
  int sum = 0;
  for (int i=0; i<length; i++)
    sum += a[i];
  return sum;
}

int findLongest(int * a, int length)
{
  int l = 0;
  int index;
  for (int i=0; i<length; i++)
  {
    if (a[i] > l)
    {
      l = a[i];
      index = i;
    }
  }
  return index;
}

int getShorter (int * a, int l, int j, int cols, int rows)
{
  int sndlong = 0;
  int min;
  for (int i=0; i<rows; i++)
  {
    int index = cols*i+j;
    if ((a[index] > sndlong) && (a[index] < l))
    {
      min = floor(log10((double)cols))+floor(log10((double)rows))+5;
      if ((a[index] < min) && (min < l))
        sndlong = min;
      else
        sndlong = a[index];
    }
  }
  if (sndlong == 0)
  {
    min = floor(log10((double)cols))+floor(log10((double)rows))+5;
    if (min < l)
      sndlong = min;
    else
      sndlong = 1;
  }
  return sndlong;
}


int * bigintmat::getwid(int maxwid)
{
  int const c = /*2**/(col-1)+1;
  if (col + c > maxwid-1) return NULL;
  int * wv = (int*)omAlloc(sizeof(int)*col*row);
  int * cwv = (int*)omAlloc(sizeof(int)*col);
  for (int j=0; j<col; j++)
  {
    cwv[j] = 0;
    for (int i=0; i<row; i++)
    {
      StringSetS("");
      n_Write(v[col*i+j], basecoeffs());
      char * tmp = StringEndS();
      const int _nl = strlen(tmp);
      wv[col*i+j] = _nl;
      if (_nl > cwv[j])
        cwv[j]=_nl;
      omFree(tmp);
    }
  }

  // Groesse verkleinern, bis < maxwid
  while (intArrSum(cwv, col)+c > maxwid)
  {
    int j = findLongest(cwv, col);
    cwv[j] = getShorter(wv, cwv[j], j, col, row);
  }
  omFree(wv);
  return cwv;
}

void bigintmat::pprint(int maxwid)
{
  if ((col==0) || (row==0))
    PrintS("");
  else
  {
    int * colwid = getwid(maxwid);
    if (colwid == NULL)
    {
      WerrorS("not enough space to print bigintmat");
      return;
    }
    char * ps;
    int slength = 0;
    for (int j=0; j<col; j++)
      slength += colwid[j]*row;
    slength += col*row+row;
    ps = (char*) omAlloc0(sizeof(char)*(slength));
    int pos = 0;
    for (int i=0; i<col*row; i++)
    {
      StringSetS("");
      n_Write(v[i], basecoeffs());
      char * ts = StringEndS();
      const int _nl = strlen(ts);
      int cj = i%col;
      if (_nl > colwid[cj])
      {
        StringSetS("");
        int ci = i/col;
        StringAppend("[%d,%d]", ci+1, cj+1);
        char * ph = StringEndS();
        int phl = strlen(ph);
        if (phl > colwid[cj])
        {
          for (int j=0; j<colwid[cj]-1; j++)
            ps[pos+j] = ' ';
          ps[pos+colwid[cj]-1] = '*';
        }
        else
        {
          for (int j=0; j<colwid[cj]-phl; j++)
            ps[pos+j] = ' ';
          for (int j=0; j<phl; j++)
            ps[pos+colwid[cj]-phl+j] = ph[j];
        }
        omFree(ph);
      }
      else  // Mit Leerzeichen auffüllen und zahl reinschreiben
      {
        for (int j=0; j<colwid[cj]-_nl; j++)
          ps[pos+j] = ' ';
        for (int j=0; j<_nl; j++)
          ps[pos+colwid[cj]-_nl+j] = ts[j];
      }
      // ", " und (evtl) "\n" einfügen
      if ((i+1)%col == 0)
      {
        if (i != col*row-1)
        {
          ps[pos+colwid[cj]] = ',';
          ps[pos+colwid[cj]+1] = '\n';
          pos += colwid[cj]+2;
        }
      }
      else
      {
        ps[pos+colwid[cj]] = ',';
        pos += colwid[cj]+1;
      }

      omFree(ts);  // Hier ts zerstören
    }
    PrintS(ps);
   // omFree(ps);
  }
}


void bigintmat::swap(int i, int j) {
  if ((i <= col) && (j <= col) && (i>0) && (j>0)) {
    number tmp;
    number t;
    for (int k=1; k<=row; k++) {
      tmp = get(k, i);
      t = get(k, j);
      set(k, i, t);
      set(k, j, tmp);
      n_Delete(&t, basecoeffs());
      n_Delete(&tmp, basecoeffs());
    }
  }
  else
    Werror("Error in swap");
}

void bigintmat::swaprow(int i, int j) {
  if ((i <= row) && (j <= row) && (i>0) && (j>0)) {
    number tmp;
    number t;
    for (int k=1; k<=col; k++) {
      tmp = get(i, k);
      t = get(j, k);
      set(i, k, t);
      set(j, k, tmp);
      n_Delete(&tmp, basecoeffs());
      n_Delete(&t, basecoeffs());
    }
  }
  else
    Werror("Error in swaprow");
}

int bigintmat::findnonzero(int i)
{
  // Setze t2 := 0. Vergleiche jeden Eintrag mit t2. Falls ungleich, gebe Index zurück, sonst teste nächsten.
  // Wurde der letzte Eintrag getestet und es war immer gleich, gebe 0 zurück.
  number t1;
  number t2 = n_Init(0, basecoeffs());
  for (int j=1; j<=col; j++) {
    t1 = get(i, j);
    if (!n_Equal(t1, t2, basecoeffs()))
    {
      n_Delete(&t1, basecoeffs());
      n_Delete(&t2, basecoeffs());
      return j;
    }
    n_Delete(&t1, basecoeffs());
  }
  n_Delete(&t2, basecoeffs());
  return 0;
}

int bigintmat::findcolnonzero(int j)
{
  // Setze t2 := 0. Vergleiche jeden Eintrag mit t2. Falls ungleich, gebe Index zurück, sonst teste nächsten.
  // Wurde der letzte Eintrag getestet und es war immer gleich, gebe 0 zurück.
  number t1;
  number t2 = n_Init(0, basecoeffs());
  for (int i=row; i>=1; i--) {
    t1 = get(i, j);
    if (!n_Equal(t1, t2, basecoeffs()))
    {
      n_Delete(&t1, basecoeffs());
      n_Delete(&t2, basecoeffs());
      return i;
    }
    n_Delete(&t1, basecoeffs());
  }
  n_Delete(&t2, basecoeffs());
  return 0;
}

void bigintmat::getcol(int j, bigintmat *a)
{
  if ((j>col) || (j<1)) {
    Werror("Error in getcol: Index out of range!");
    return;
  }
  if (((a->rows() != row) || (a->cols() != 1)) && ((a->rows() != 1) || (a->cols() != row))) {
    Werror("Error in getcol. Dimensions must agree!");
    return;
  }
  if (!nCoeffs_are_equal(basecoeffs(), a->basecoeffs())) {
    nMapFunc f = n_SetMap(basecoeffs(), a->basecoeffs());
    number t1, t2;
    for (int i=1; i<=row;i++) {
    t1 = get(i,j);
    t2 = f(t1, basecoeffs(), a->basecoeffs());
    a->set(i-1,t1);
    n_Delete(&t1, basecoeffs());
    n_Delete(&t2, a->basecoeffs());
  }
    return;
  }
  number t1;
  for (int i=1; i<=row;i++) {
    t1 = get(i,j);
    a->set(i-1,t1);
    n_Delete(&t1, basecoeffs());
  }
}

void bigintmat::getrow(int i, bigintmat *a)
{
  if ((i>row) || (i<1)) {
    Werror("Error in getrow: Index out of range!");
    return;
  }
  if (((a->rows() != 1) || (a->cols() != col)) && ((a->rows() != col) || (a->cols() != 1))) {
    Werror("Error in getrow. Dimensions must agree!");
    return;
  }
  if (!nCoeffs_are_equal(basecoeffs(), a->basecoeffs())) {
    nMapFunc f = n_SetMap(basecoeffs(), a->basecoeffs());
    number t1, t2;
    for (int j=1; j<=col;j++) {
      t1 = get(i,j);
      t2 = f(t1, basecoeffs(), a->basecoeffs());
      a->set(j-1,t2);
      n_Delete(&t1, basecoeffs());
      n_Delete(&t2, a->basecoeffs());
    }
    return;
  } 
  number t1;
  for (int j=1; j<=col;j++) {
    t1 = get(i,j);
    a->set(j-1,t1);
    n_Delete(&t1, basecoeffs());
  }
}


void bigintmat::setcol(int j, bigintmat *m)
{
   if ((j>col) || (j<1)) {
    Werror("Error in setcol: Index out of range!");
    return;
  }
  if (((m->rows() != row) || (m->cols() != 1)) && ((m->rows() != 1) || (m->cols() != row))) {
    Werror("Error in setcol. Dimensions must agree!");
    return;
  }
  if (!nCoeffs_are_equal(basecoeffs(), m->basecoeffs())) {
    nMapFunc f = n_SetMap(m->basecoeffs(), basecoeffs());
    number t1,t2;
    for (int i=1; i<=row; i++) {
      t1 = m->get(i-1);
      t2 = f(t1, m->basecoeffs(),basecoeffs());
      set(i, j, t2);
      n_Delete(&t2, basecoeffs());
      n_Delete(&t1, m->basecoeffs());
    }
    return;
  }
  number t1;
  for (int i=1; i<=row; i++) {
      t1 = m->get(i-1);
      set(i, j, t1);
      n_Delete(&t1, basecoeffs());
  }
}

void bigintmat::setrow(int j, bigintmat *m) {
  if ((j>row) || (j<1)) {
    Werror("Error in setrow: Index out of range!");
    return;
  }
  if (((m->rows() != 1) || (m->cols() != col)) && ((m->rows() != col) || (m->cols() != 1))) {
    Werror("Error in setrow. Dimensions must agree!");
    return;
  }
  if (!nCoeffs_are_equal(basecoeffs(), m->basecoeffs())) {
    nMapFunc f = n_SetMap(m->basecoeffs(), basecoeffs());
    number tmp1,tmp2;
  for (int i=1; i<=col; i++) {
      tmp1 = m->get(i-1);
      tmp2 = f(tmp1, m->basecoeffs(),basecoeffs());
      set(j, i, tmp2);
      n_Delete(&tmp2, basecoeffs());
      n_Delete(&tmp1, m->basecoeffs());
  }
    return;
  }
  number tmp;
  for (int i=1; i<=col; i++) {
      tmp = m->get(i-1);
      set(j, i, tmp);
      n_Delete(&tmp, basecoeffs());
  }
  
}

bool bigintmat::add(bigintmat *b)
{
  if ((b->rows() != row) || (b->cols() != col)) {
    Werror("Error in bigintmat::add. Dimensions do not agree!");
    return false;
  }
  if (!nCoeffs_are_equal(basecoeffs(), b->basecoeffs())) {
    Werror("Error in bigintmat::add. coeffs do not agree!");
    return false;
  }
  bigintmat *sum = bimAdd(this, b);
  number t1;
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      t1 = sum->get(i, j);
      set(i, j, t1);
      n_Delete(&t1, basecoeffs());
    }
  }
  delete sum;
  return true;
}

bool bigintmat::sub(bigintmat *b)
{
 if ((b->rows() != row) || (b->cols() != col)) {
    Werror("Error in bigintmat::sub. Dimensions do not agree!");
    return false;
  }
  if (!nCoeffs_are_equal(basecoeffs(), b->basecoeffs())) {
    Werror("Error in bigintmat::sub. coeffs do not agree!");
    return false;
  }
  bigintmat *dif = bimSub(this, b);
  number t1;
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      t1 = dif->get(i, j);
      set(i, j, t1);
      n_Delete(&t1, basecoeffs());
    }
  }
  delete dif;
  return true;
}

bool bigintmat::skalmult(number b, coeffs c)
{
  if (!nCoeffs_are_equal(c, basecoeffs()))
  {
    Werror("Wrong coeffs\n");
    return false;
  }
  number t1, t2;
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      t1 = get(i, j);
      t2 = n_Mult(t1, b, basecoeffs());
      set(i, j, t2);
      n_Delete(&t1, basecoeffs());
      n_Delete(&t2, basecoeffs());
    }
  }
  return true;
}

bool bigintmat::addcol(int i, int j, number a, coeffs c)
{
  if ((i>col) || (j>col) || (i<1) || (j<1)) {
    Werror("Error in addcol: Index out of range!");
    return false;
  }
  if (!nCoeffs_are_equal(c, basecoeffs())) {
    Werror("Error in addcol: coeffs do not agree!");
    return false;
  }
  number t1, t2, t3, t4;
  for (int k=1; k<=row; k++)
  {
    t1 = get(k, j);
    t2 = get(k, i);
    t3 = n_Mult(t1, a, basecoeffs());
    t4 = n_Add(t3, t2, basecoeffs());
    set(k, i, t4);
    n_Delete(&t3, basecoeffs());
    n_Delete(&t4, basecoeffs());
    n_Delete(&t1, basecoeffs());
    n_Delete(&t2, basecoeffs());
  }
  return true;
}

bool bigintmat::addrow(int i, int j, number a, coeffs c)
{
  if ((i>row) || (j>row) || (i<1) || (j<1)) {
    Werror("Error in addrow: Index out of range!");
    return false;
  }
  if (!nCoeffs_are_equal(c, basecoeffs())) {
    Werror("Error in addrow: coeffs do not agree!");
    return false;
  }
  number t1, t2, t3, t4;
  for (int k=1; k<=col; k++)
  {
    t1 = get(j, k);
    t2 = get(i, k);
    t3 = n_Mult(t1, a, basecoeffs());
    t4 = n_Add(t3, t2, basecoeffs());
    set(i, k, t4);
    n_Delete(&t1, basecoeffs());
    n_Delete(&t2, basecoeffs());
    n_Delete(&t3, basecoeffs());
    n_Delete(&t4, basecoeffs());
  }
  return true;
}

void bigintmat::colskalmult(int i, number a, coeffs c) {
  if ((i>=1) && (i<=col) && (nCoeffs_are_equal(c, basecoeffs()))) {
    number t, tmult;
    for (int j=1; j<=row; j++) {
      t = get(j, i);
      tmult = n_Mult(a, t, basecoeffs());
      set(j, i, tmult);
      n_Delete(&t, basecoeffs());
      n_Delete(&tmult, basecoeffs());
    }
  }
  else
    Werror("Error in colskalmult");
}

void bigintmat::rowskalmult(int i, number a, coeffs c) {
  if ((i>=1) && (i<=row) && (nCoeffs_are_equal(c, basecoeffs()))) {
    number t, tmult;
    for (int j=1; j<=col; j++) {
      t = get(i, j);
      tmult = n_Mult(a, t, basecoeffs());
      set(i, j, tmult);
      n_Delete(&t, basecoeffs());
      n_Delete(&tmult, basecoeffs());
    }
  }
  else
    Werror("Error in rowskalmult");
}

void bigintmat::concatrow(bigintmat *a, bigintmat *b) {
  int ay = a->cols();
  int ax = a->rows();
  int by = b->cols();
  int bx = b->rows();
  number tmp;
  if (!((col == ay) && (col == by) && (ax+bx == row))) {
    Werror("Error in concatrow. Dimensions must agree!");
    return;
  }
  if (!(nCoeffs_are_equal(a->basecoeffs(), basecoeffs()) && nCoeffs_are_equal(b->basecoeffs(), basecoeffs()))) {
    Werror("Error in concatrow. coeffs do not agree!");
    return;
  }
  for (int i=1; i<=ax; i++) {
    for (int j=1; j<=ay; j++) {
      tmp = a->get(i,j);
      set(i, j, tmp);
      n_Delete(&tmp, basecoeffs());
    }
  }
  for (int i=1; i<=bx; i++) {
    for (int j=1; j<=by; j++) {
      tmp = b->get(i,j);
      set(i+ax, j, tmp);
      n_Delete(&tmp, basecoeffs());
    }
  }
} 

void bigintmat::concatcol (bigintmat *a, bigintmat *b) {
  int ay = a->cols();
  int ax = a->rows();
  int by = b->cols();
  int bx = b->rows();
  number tmp;
  if (!((row == ax) && (row == bx) && (ay+by==col))) {
    Werror("Error in concatcol");
    return;
  }
  if (!(nCoeffs_are_equal(a->basecoeffs(), basecoeffs()) && nCoeffs_are_equal(b->basecoeffs(), basecoeffs()))) {
    Werror("Error in concatcol. coeffs do not agree!");
    return;
  }
  for (int i=1; i<=ax; i++) {
    for (int j=1; j<=ay; j++) {
      tmp = a->get(i,j);
      set(i, j, tmp);
      n_Delete(&tmp, basecoeffs());
    }
  }
  for (int i=1; i<=bx; i++) {
    for (int j=1; j<=by; j++) {
      tmp = b->get(i,j);
      set(i, j+ay, tmp);
      n_Delete(&tmp, basecoeffs());
    }
  }
}

void bigintmat::splitrow(bigintmat *a, bigintmat *b) {
  int ay = a->cols();
  int ax = a->rows();
  int by = b->cols();
  int bx = b->rows();
	number tmp;
  if (!(ax + bx == row)) {
    Werror("Error in splitrow. Dimensions must agree!");
  }
  else if (!((col == ay) && (col == by))) {
    Werror("Error in splitrow. Dimensions must agree!");
  }
  else if (!(nCoeffs_are_equal(a->basecoeffs(), basecoeffs()) && nCoeffs_are_equal(b->basecoeffs(), basecoeffs()))) {
    Werror("Error in splitrow. coeffs do not agree!");
  }
  else {
    for(int i = 1; i<=ax; i++) {
      for(int j = 1; j<=ay;j++) {
        tmp = get(i,j);
        a->set(i,j,tmp);
        n_Delete(&tmp, basecoeffs());
      }
    }
    for (int i =1; i<=bx; i++) {
      for (int j=1;j<=col;j++) {
        tmp = get(i+ax, j);
        b->set(i,j,tmp);
        n_Delete(&tmp, basecoeffs());
      }
    }
  }
}

void bigintmat::splitcol(bigintmat *a, bigintmat *b) {
  int ay = a->cols();
  int ax = a->rows();
  int by = b->cols();
  int bx = b->rows();
  number tmp;
  if (!((row == ax) && (row == bx))) {
    Werror("Error in splitcol. Dimensions must agree!");
  }
  else if (!(ay+by == col)) {
    Werror("Error in splitcol. Dimensions must agree!");
  }
  else if (!(nCoeffs_are_equal(a->basecoeffs(), basecoeffs()) && nCoeffs_are_equal(b->basecoeffs(), basecoeffs()))) {
    Werror("Error in splitcol. coeffs do not agree!");
  }
  else {
    for (int i=1; i<=ax; i++) {
      for (int j=1; j<=ay; j++) {
        tmp = get(i,j);
        a->set(i,j,tmp);
        n_Delete(&tmp, basecoeffs());
      }
    }
    for (int i=1; i<=bx; i++) {
      for (int j=1; j<=by; j++) {
        tmp = get(i,j+ay);
        b->set(i,j,tmp);
        n_Delete(&tmp, basecoeffs());
      }
    }
  }
}

void bigintmat::splitcol(bigintmat *a, int i) {
  number tmp;
  if ((a->rows() != row) || (a->cols()+i-1 > col) || (i<1)) {
    Werror("Error in splitcol. Dimensions must agree!");
    return;
  }
  if (!(nCoeffs_are_equal(a->basecoeffs(), basecoeffs()))) {
    Werror("Error in splitcol. coeffs do not agree!");
    return;
  }
  int width = a->cols();
  for (int j=1; j<=width; j++) {
    for (int k=1; k<=row; k++) {
      tmp = get(k, j+i-1);
      a->set(k, j, tmp);
      n_Delete(&tmp, basecoeffs());
    }
  }
}

void bigintmat::splitrow(bigintmat *a, int i) {
  number tmp;
  if ((a->cols() != col) || (a->rows()+i-1 > row) || (i<1)) {
    Werror("Error in Marco-splitrow");
    return;
  }
  
  if (!(nCoeffs_are_equal(a->basecoeffs(), basecoeffs()))) {
    Werror("Error in splitrow. coeffs do not agree!");
    return;
  }
  int height = a->rows();
  for (int j=1; j<=height; j++) {
    for (int k=1; k<=col; k++) {
      tmp = view(j+i-1, k);
      a->set(j, k, tmp);
    }
  }
}

bool bigintmat::copy(bigintmat *b)
{
  if ((b->rows() != row) || (b->cols() != col)) {
    Werror("Error in bigintmat::copy. Dimensions do not agree!");
    return false;
  }
  if (!nCoeffs_are_equal(basecoeffs(), b->basecoeffs())) {
    Werror("Error in bigintmat::copy. coeffs do not agree!");
    return false;
  }
  number t1;
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      t1 = b->view(i, j);
      set(i, j, t1);
    }
  }
  return true;
}

int bigintmat::isOne() {
  coeffs r = basecoeffs();
  if (row==col) {
    for (int i=1; i<=row; i++) {
      for (int j=1; j<=col; j++) {
        if (i==j) {
          if (!n_IsOne(view(i, j), r))
            return 0;
        }
        else {
          if (!n_IsZero(view(i,j), r))
            return 0;
        }
      }
    }
  }
  return 1;
}


void bigintmat::one() {
  if (row==col) {
    number one = n_Init(1, basecoeffs()),
           zero = n_Init(0, basecoeffs());
    for (int i=1; i<=row; i++) {
      for (int j=1; j<=col; j++) {
        if (i==j) {
          set(i, j, one);
        }
        else {
          set(i, j, zero);
        }
      }
    }
    n_Delete(&one, basecoeffs());
    n_Delete(&zero, basecoeffs());
  }
}

void bigintmat::zero() {
  number tmp = n_Init(0, basecoeffs());
  for (int i=1; i<=row; i++) {
    for (int j=1; j<=col; j++) {
      set(i, j, tmp);
    }
  }
  n_Delete(&tmp,basecoeffs());
}


//used in the det function. No idea what it does.
//looks like it return the submatrix where the i-th row
//and j-th column has been removed in the LaPlace generic
//determinant algorithm
bigintmat *bigintmat::elim(int i, int j)
{
  if ((i<=0) || (i>row) || (j<=0) || (j>col))
    return NULL;
  int cx, cy;
  cx=1;
  cy=1;
  number t;
  bigintmat *b = new bigintmat(row-1, col-1, basecoeffs());
  for (int k=1; k<=row; k++) {
    if (k!=i)
    {
      cy=1;
      for (int l=1; l<=col; l++)
      {
        if (l!=j)
        {
          t = get(k, l);
          b->set(cx, cy, t);
          n_Delete(&t, basecoeffs());
          cy++;
        }
      }
      cx++;
    }
  }
  return b;
}
#endif


//returns d such that a/d is the inverse of the input
//TODO: make work for p not prime using the euc stuff.
//long term: rewrite for Z using p-adic lifting
//and Dixon. Possibly even the sparse recent Storjohann stuff
number bigintmat::pseudoinv(bigintmat *a) {

  // Falls Matrix über reellen Zahlen nicht invertierbar, breche ab
 assume((a->rows() == row) && (a->rows() == a->cols()) && (row == col));
    
 number det = this->det(); //computes the HNF, so should e reused.
 if ((n_IsZero(det, basecoeffs())))
    return det;

 // Hänge Einheitsmatrix über Matrix und wendet HNF an. An Stelle der Einheitsmatrix steht im Ergebnis die Transormationsmatrix dazu
  a->one();
  bigintmat *m = new bigintmat(2*row, col, basecoeffs());
  m->concatrow(a,this);
  m->hnf();
  // Arbeite weiterhin mit der zusammengehängten Matrix
  // Laufe durch die Diagonalelemente, und multipliziere jede Spalte rechts davon damit, speichere aber den alten Eintrag der Spalte, temp, der in der Zeile des Diagonalelements liegt, zwischen. Dann addiere das -temp-Fache der Diagonalspalte zur entsprechenenden Spalte rechts davon. Dadurch entsteht überall rechts der Diagonalen eine 0
  number diag;
  number temp, ttemp;
  for (int i=1; i<=col; i++) {
    diag = m->get(row+i, i);
    for (int j=i+1; j<=col; j++) {
      temp = m->get(row+i, j);
      m->colskalmult(j, diag, basecoeffs());
      temp = n_Neg(temp, basecoeffs());
      m->addcol(j, i, temp, basecoeffs());
      n_Delete(&temp, basecoeffs());
    }
    n_Delete(&diag, basecoeffs());
  }
  // Falls wir nicht modulo n arbeiten, können wir die Spalten durch den ggT teilen, um die Einträge kleiner zu bekommen
  // Bei Z/n sparen wir uns das, da es hier sinnlos ist
  number g;
  number gcd;
  for (int j=1; j<=col; j++) {
    g = n_Init(0, basecoeffs());
    for (int i=1; i<=2*row; i++) {
      temp = m->get(i,j);
      gcd = n_Gcd(g, temp, basecoeffs());
      n_Delete(&g, basecoeffs());
      n_Delete(&temp, basecoeffs());
      g = n_Copy(gcd, basecoeffs());
      n_Delete(&gcd, basecoeffs());
    }
    if (!(n_IsOne(g, basecoeffs())))
      m->colskaldiv(j, g);
    n_Delete(&g, basecoeffs());
  }
  
  // Nun müssen die Diagonalelemente durch Spaltenmultiplikation gleich gesett werden. Bei Z können wir mit dem kgV arbeiten, bei Z/n bringen wir jedes Diagonalelement auf 1 (wir arbeiten immer mit n = Primzahl. Für n != Primzahl muss noch an anderen Stellen etwas geändert werden)
  
  g = n_Init(0, basecoeffs());
  number prod = n_Init(1, basecoeffs());
  for (int i=1; i<=col; i++) {
    gcd = n_Gcd(g, m->get(row+i, i), basecoeffs());
    n_Delete(&g, basecoeffs());
    g = n_Copy(gcd, basecoeffs());
    n_Delete(&gcd, basecoeffs());
    ttemp = n_Copy(prod, basecoeffs());
    temp = m->get(row+i, i);
    n_Delete(&prod, basecoeffs());
    prod = n_Mult(ttemp, temp, basecoeffs());
    n_Delete(&ttemp, basecoeffs());
    n_Delete(&temp, basecoeffs());
  }
  number lcm;
  lcm = n_Div(prod, g, basecoeffs());
  for (int j=1; j<=col; j++) {
    ttemp = m->get(row+j,j);
    temp = n_QuotRem(lcm, ttemp, NULL, basecoeffs());
    m->colskalmult(j, temp, basecoeffs());
    n_Delete(&ttemp, basecoeffs());
    n_Delete(&temp, basecoeffs());
  }
  n_Delete(&lcm, basecoeffs());
  n_Delete(&prod, basecoeffs());
  
  number divisor = m->get(row+1, 1);
  m->splitrow(a, 1);
  delete m;
  n_Delete(&det, basecoeffs());
  return divisor;
}

number bigintmat::trace()
{
  assume (col == row);
  number t = get(1,1),
         h;
  coeffs r = basecoeffs();
  for(int i=2; i<= col; i++) {
    h = n_Add(t, view(i,i), r);
    n_Delete(&t, r);
    t = h;
  }
  return t;
}

number bigintmat::det()
{
  assume (row==col);

  if (col == 1)
    return get(1, 1);
  // should work as well in Z/pZ of type n_Zp?
  // elies on XExtGcd and the other euc. functinos.
  if ( getCoeffType(basecoeffs())== n_Z || getCoeffType(basecoeffs() )== n_Zn) {
    return hnfdet();
  }
  number sum = n_Init(0, basecoeffs());
  number t1, t2, t3, t4;
  bigintmat *b;
  for (int i=1; i<=row; i++) {
    b = elim(i, 1);
    t1 = get(i, 1);
    t2 = b->det();
    t3 = n_Mult(t1, t2, basecoeffs());
    t4 = n_Copy(sum, basecoeffs());
    n_Delete(&sum, basecoeffs());
    if ((i+1)>>1<<1==(i+1))
      sum = n_Add(t4, t3, basecoeffs());
    else
      sum = n_Sub(t4, t3, basecoeffs());
    n_Delete(&t1, basecoeffs());
    n_Delete(&t2, basecoeffs());
    n_Delete(&t3, basecoeffs());
    n_Delete(&t4, basecoeffs());
  }
  return sum;
}

number bigintmat::hnfdet()
{
  assume (col == row);

  if (col == 1)
    return get(1, 1);
  bigintmat *m = new bigintmat(this);
  m->hnf();
  number prod = n_Init(1, basecoeffs());
  number temp, temp2;
  for (int i=1; i<=col; i++) {
    temp = m->get(i, i);
    temp2 = n_Mult(temp, prod, basecoeffs());
    n_Delete(&prod, basecoeffs());
    prod = temp2;
    n_Delete(&temp, basecoeffs());
  }
  return prod;
}

#if 0
//reduce the row i modulo p, not used
void bigintmat::rowmod(int i, number p, coeffs c)
{
  if ((i>row) || (i<1)) {
    Werror("Error in addrow: Index out of range!");
    return;
  }
  coeffs coe = numbercoeffs(p, c);
  nMapFunc f1 = n_SetMap(basecoeffs(), coe);
  nMapFunc f2 = n_SetMap(coe, basecoeffs());
  number a;
  number b;
  for (int j=1; j<=col; j++) {
    a = get(i, j);
    b = f1(a, basecoeffs(), coe);
    n_Delete(&a, basecoeffs());
    a = f2(b, coe, basecoeffs());
    n_Delete(&b, coe);
    set(i, j, a);
    n_Delete(&a, basecoeffs());
  }
}
#endif


void bigintmat::hnf()
{
  // Keine 100%ige HNF, da Einträge rechts von Diagonalen auch kleiner 0 sein können
  // Laufen von unten nach oben und von links nach rechts

  if (getCoeffType(basecoeffs()) == n_Q) {
    coeffs s = nInitChar(n_Z, NULL);
    bigintmat *m = bimChangeCoeff(this, s);
    ::Print("mat over Z is \n");
    ::Print("%s\n(%d x %d)\n", m->String(), m->rows(), m->cols());
    m->hnf();
    ::Print("hnf over Z is\n%s\n", m->String());
  }

  int i = rows();
  int j = cols();
  number q = n_Init(0, basecoeffs());
  number one = n_Init(1, basecoeffs());
  number minusone = n_Init(-1, basecoeffs());
  number tmp1 = n_Init(0, basecoeffs());
  number tmp2 = n_Init(0, basecoeffs());
  number co1 = n_Init(0, basecoeffs());
  number co2 = n_Init(0, basecoeffs());
  number co3 = n_Init(0, basecoeffs());
  number co4 = n_Init(0, basecoeffs());
  number ggt = n_Init(0, basecoeffs());

  while ((i>0) && (j>0)) {
    // Falls erstes Nicht-Null-Element in Zeile i nicht existiert, oder hinter Spalte j vorkommt, gehe in nächste Zeile 
    if ((findnonzero(i)==0) || (findnonzero(i)>j)) {
      i--;
    } else {
      // Laufe von links nach rechts durch die Zeile:
      for (int l=1; l<=j-1; l++) {
        n_Delete(&tmp1, basecoeffs());
        tmp1 = get(i, l);
        // Falls Eintrag (im folgenden x genannt) gleich 0, gehe eine Spalte weiter. Ansonsten...
        if (!n_IsZero(tmp1, basecoeffs())) {
          n_Delete(&tmp2, basecoeffs());
          tmp2 = get(i, l+1);
          // Falls Eintrag (i.f. y g.) rechts daneben gleich 0, tausche beide Spalten, sonst...
          if (!n_IsZero(tmp2, basecoeffs())) {
            n_Delete(&ggt, basecoeffs());
            ggt = n_XExtGcd(tmp1, tmp2, &co1, &co2, &co3, &co4, basecoeffs());
            // Falls x=ggT(x, y), tausche die beiden Spalten und ziehe die (neue) rechte Spalte so häufig von der linken ab, dass an der ehemaligen Stelle von x nun eine 0 steht. Dazu:
            if (n_Equal(tmp1, ggt, basecoeffs())) { 
              swap(l, l+1);
              n_Delete(&q, basecoeffs());
              q = n_Div(tmp2, ggt, basecoeffs());
              q = n_Neg(q, basecoeffs());
              // Dann addiere das -q-fache der (neuen) rechten Spalte zur linken dazu. Damit erhalten wir die gewünschte 0

              addcol(l, l+1, q, basecoeffs());
              n_Delete(&q, basecoeffs());
              
            }
            else if (n_Equal(tmp1, minusone, basecoeffs())) {
              // Falls x=-1, so ist x=-ggt(x, y). Dann gehe wie oben vor, multipliziere aber zuerst die neue rechte Spalte (die mit x) mit -1
              // Die Berechnung von q (=y/ggt) entfällt, da ggt=1
              swap(l, l+1);
              colskalmult(l+1, minusone, basecoeffs());
              tmp2 = n_Neg(tmp2, basecoeffs());
              addcol(l, l+1, tmp2, basecoeffs());
            }
            else {
              // CF: use the 2x2 matrix (co1, co2)(co3, co4) to
              // get the gcd in position and the 0 in the other:
#ifdef CF_DEB             
              Print("applying trafo\n");
              StringSetS("");
              n_Write(co1, basecoeffs()); StringAppendS("\t");
              n_Write(co2, basecoeffs()); StringAppendS("\t");
              n_Write(co3, basecoeffs()); StringAppendS("\t");
              n_Write(co4, basecoeffs()); StringAppendS("\t");
              Print("%s\nfor l=%d\n", StringEndS(), l);
              Print("to %s\n", String());
#endif
              coltransform(l, l+1, co3, co4, co1, co2);
#ifdef CF_DEB              
              Print("gives  %s\n", String());
#endif
            }
            n_Delete(&co1, basecoeffs());
            n_Delete(&co2, basecoeffs());
            n_Delete(&co3, basecoeffs());
            n_Delete(&co4, basecoeffs());
          }
          else {
            swap(l, l+1);
          }
          // Dann betrachte die vormals rechte Spalte als neue linke, und die rechts daneben als neue rechte.
        }
      }
      
      tmp1 = get(i, j);
      // normalize by units:
      if (!n_IsZero(tmp1, basecoeffs())) {
        number u = n_GetUnit(tmp1, basecoeffs());
        if (!n_IsOne(u, basecoeffs())) {
          colskalmult(j, u, basecoeffs());
        }
        n_Delete(&u, basecoeffs());
      }
      // Zum Schluss mache alle Einträge rechts vom Diagonalelement betragsmäßig kleiner als dieses
      // Durch Änderung könnte man auch erreichen, dass diese auch nicht negativ sind
      for (int l=j+1; l<=col; l++) {
        tmp1 = get(i, l);
        tmp2 = get(i, j);
        n_Delete(&q, basecoeffs());
        q = n_QuotRem(tmp1, tmp2, NULL, basecoeffs());
        q = n_Neg(q, basecoeffs());
        addcol(l, j, q, basecoeffs());
      }
      i--;
      j--;
      // Dann betrachte die Zeile darüber und gehe dort wie vorher vor
    }
  }
  n_Delete(&q, basecoeffs());
  n_Delete(&tmp1, basecoeffs());
  n_Delete(&tmp2, basecoeffs());
  n_Delete(&ggt, basecoeffs());
  n_Delete(&co1, basecoeffs());
  n_Delete(&co2, basecoeffs());
  n_Delete(&one, basecoeffs());
}

void bigintmat::modhnf2(number p, coeffs c)
{
  //CF: seems to do the same as above - until the last step, when the
  //  entire matrix seems to e "reduced mod p".
  //
  // Keine 100%ige HNF, da Einträge rechts von Diagonalen auch kleiner 0 sein können
  // Laufen von unten nach oben und von links nach rechts
  int i = rows();
  int j = cols();
  number invggt;
  number q;
  number minusone = n_Init(-1, basecoeffs());
  number tmp1 = n_Init(0, basecoeffs());
  number tmp2 = n_Init(0, basecoeffs());
  number co1 = n_Init(0, basecoeffs());
  number co2 = n_Init(0, basecoeffs());
  number ggt = n_Init(0, basecoeffs());

  while ((i>0) && (j>0)) {
    // Falls erstes Nicht-Null-Element in Zeile i nicht existiert, oder hinter Spalte j vorkommt, gehe in nächste Zeile 
    if ((findnonzero(i)==0) || (findnonzero(i)>j))
      i--;
    else {
      // Laufe von links nach rechts durch die Zeile:
      for (int l=1; l<=j-1; l++) {
        tmp1 = get(i, l);
        // Falls Eintrag (im folgenden x genannt) gleich 0, gehe eine Spalte weiter. Ansonsten...
        if (!n_IsZero(tmp1, basecoeffs())) {
          tmp2 = get(i, l+1);
          // Falls Eintrag (i.f. y g.) rechts daneben gleich 0, tausche beide Spalten, sonst...
          if (!n_IsZero(tmp2, basecoeffs())) {
            ggt = n_ExtGcd(tmp1, tmp2, &co1, &co2, basecoeffs());
            // Falls x=ggT(x, y), tausche die beiden Spalten und ziehe die (neue) rechte Spalte so häufig von der linken ab, dass an der ehemaligen Stelle von x nun eine 0 steht. Dazu:
            if (n_Equal(tmp1, ggt, basecoeffs())) { 
              swap(l, l+1);
              // Falls wir modulo p rechnen: Multipliziere y mit dem Inversen von ggt, und speichere als q zwischen
              if (getCoeffType(basecoeffs())==n_Zn) {
                invggt = n_Invers(ggt, basecoeffs());
                q = n_Mult(tmp2, invggt, basecoeffs());
                n_Delete(&invggt, basecoeffs());
              }
              else // Ansonsten dividieren wir y durch ggt (ganzzahldivision, geht aber auf)
                q = n_IntDiv(tmp2, ggt, basecoeffs());
              q = n_Neg(q, basecoeffs());
              // Dann addiere das -q-fache der (neuen) rechten Spalte zur linken dazu. Damit erhalten wir die gewünschte 0
              addcol(l, l+1, q, basecoeffs());
              n_Delete(&q, basecoeffs());
            }
            else if (n_Equal(tmp1, minusone, basecoeffs())) {
              // Falls x=-1, so ist x=-ggt(x, y). Dann gehe wie oben vor, multipliziere aber zuerst die neue rechte Spalte (die mit x) mit -1
              // Die Berechnung von q (=y/ggt) entfällt, da ggt=1
              swap(l, l+1);
              colskalmult(l+1, minusone, basecoeffs());
              tmp2 = n_Neg(tmp2, basecoeffs());
              addcol(l, l+1, tmp2, basecoeffs());
            }
            else {
              // Sonst bringe durch Spaltenoperationen und Koeffizienten des erweiterten ggT den Eintrag an Stelle von y auf den ggt, und ziehe diese Spalte dann entsprechend häufig von der links daneben liegenden ab, um an Stelle des x eine 0 zu erhalten
              colskalmult(l+1, co2, basecoeffs());
              addcol(l+1, l, co1, basecoeffs());
              if (getCoeffType(basecoeffs())==n_Zn) {
                invggt = n_Invers(ggt, basecoeffs());
                q = n_Mult(tmp1, invggt, basecoeffs());
                n_Delete(&invggt, basecoeffs());
              }
              else {
                q = n_IntDiv(tmp1, ggt, basecoeffs());
              }
              q = n_Neg(q, basecoeffs());
              addcol(l, l+1, q, basecoeffs());
              n_Delete(&q, basecoeffs());
            }
          }
          else {
            swap(l, l+1);
          }
          // Dann betrachte die vormals rechte Spalte als neue linke, und die rechts daneben als neue rechte.
        }
      }
      n_Delete(&tmp1, basecoeffs());
      n_Delete(&tmp2, basecoeffs());
      n_Delete(&ggt, basecoeffs());
      n_Delete(&co1, basecoeffs());
      n_Delete(&co2, basecoeffs());
      tmp1 = get(i, j);
      // Falls wir in Z arbeiten, mache alle Diagonalelemente durch Spaltenmultiplikation mit -1 positiv
      if (getCoeffType(basecoeffs()) != n_Zn) {
        if (!(n_GreaterZero(tmp1, basecoeffs()) || n_IsZero(tmp1, basecoeffs()))) {
          tmp2 = n_Init(-1, basecoeffs());
          colskalmult(j, tmp2, basecoeffs());
          n_Delete(&tmp2, basecoeffs());
        }
      }
      // Zum Schluss mache alle Einträge rechts vom Diagonalelement betragsmäßig kleiner als dieses
      // Durch Änderung könnte man auch erreichen, dass diese auch nicht negativ sind
      for (int l=j+1; l<=col; l++) {
        tmp1 = get(i, l);
        tmp2 = get(i, j);
        if (getCoeffType(basecoeffs())==n_Zn) {//CF crap
          invggt = n_Invers(tmp2, basecoeffs());
          q = n_Mult(tmp1, invggt, basecoeffs());
          n_Delete(&invggt, basecoeffs());
        }
        else
          q = n_QuotRem(tmp1, tmp2, NULL, basecoeffs());
        q = n_Neg(q, basecoeffs());
        addcol(l, j, q, basecoeffs());
        n_Delete(&q, basecoeffs());
      }
      mod(p, c);
      i--;
      j--;
      // Dann betrachte die Zeile darüber und gehe dort wie vorher vor
    }
  }
}

bigintmat * bimChangeCoeff(bigintmat *a, coeffs cnew)
{
  coeffs cold = a->basecoeffs();
  bigintmat *b = new bigintmat(a->rows(), a->cols(), cnew);
  // Erzeugt Karte von alten coeffs nach neuen
  nMapFunc f = n_SetMap(cold, cnew);
  number t1;
  number t2;
  // Wende Karte auf jeden Eintrag von a an und schreibe Ergebnis in b
  for (int i=1; i<=a->rows(); i++)
  {
    for (int j=1; j<=a->cols(); j++)
    {
      t1 = a->get(i, j);
      t2 = f(t1, cold, cnew);
      b->set(i, j, t2);
      n_Delete(&t1, cold);
      n_Delete(&t2, cnew);
    }
  }
  return b;
}

//a lift of the HNF in Z/pZ
bigintmat * bigintmat::modhnf(number p, coeffs c)
{
  coeffs coe = numbercoeffs(p, c);
  bigintmat *m = bimChangeCoeff(this, coe);
  m->hnf();
  bigintmat *a = bimChangeCoeff(m, basecoeffs());
  delete m;
  return a;
}


//for p prime: echelon form (gaussian elemination) in Z/pZ
//should be changed to work over a field and do the coercion
//elsewhere (in the calling function)
bigintmat *bigintmat::modgauss(number p, coeffs c) {
  bigintmat *m = this->inpmod(p, c);
  number temp, temp2;
  number inv;
  int inttemp;
  for (int i=1; i<=row; i++) {
    temp = m->get(i,i);
    inttemp = m->findcolnonzero(i);
    if (n_IsZero(temp,m->basecoeffs()) && (inttemp > i)) 
      m->swaprow(i, inttemp);
    n_Delete(&temp,m->basecoeffs());
    temp = m->get(i,i);
    if (!n_IsZero(temp,m->basecoeffs())) {
      if (!n_IsOne(temp,m->basecoeffs())) {
        inv = n_Invers(temp, m->basecoeffs());
        m->rowskalmult(i, inv, m->basecoeffs());
        n_Delete(&inv, m->basecoeffs());
      }
      for (int j=i+1; j<=row; j++) {
        temp2 = m->get(j,i);
        temp2 = n_Neg(temp2, m->basecoeffs());
        m->addrow(j, i, temp2, m->basecoeffs());
        n_Delete(&temp2, m->basecoeffs());
      }
    }
    n_Delete(&temp, m->basecoeffs());
  }
  // Reduzieren der Einträge über den Diagonalen
  for (int i=row; i>=1; i--) {
    temp = m->get(i,i);
    if (!n_IsZero(temp, m->basecoeffs())) {
      for (int j=i-1; j>=1; j--) {
        temp2 = m->get(j,i);
        temp2 = n_Neg(temp2, m->basecoeffs());
        m->addrow(j, i, temp2, m->basecoeffs());
        n_Delete(&temp2, m->basecoeffs());
      }
    }
    n_Delete(&temp, m->basecoeffs());
  }
  return bimChangeCoeff(m, basecoeffs());
}

//exactly divide matrix by b
void bigintmat::skaldiv(number b)
{
  number tmp1, tmp2;
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      tmp1 = view(i, j);
      tmp2 = n_Div(tmp1, b, basecoeffs());
      rawset(i, j, tmp2);
    }
  }
}

//exactly divide col j by b
void bigintmat::colskaldiv(int j, number b)
{
  number tmp1, tmp2;
  for (int i=1; i<=row; i++)
  {
    tmp1 = view(i, j);
    tmp2 = n_Div(tmp1, b, basecoeffs());
    rawset(i, j, tmp2);
  }
}

// col(j, k) <- col(j,k)*matrix((a, c)(b, d))
// mostly used internally in the hnf and Howell stuff
void bigintmat::coltransform(int j, int k, number a, number b, number c, number d)
{
  number tmp1, tmp2, tmp3, tmp4;
  for (int i=1; i<=row; i++)
  {
    tmp1 = get(i, j);
    tmp2 = get(i, k);
    tmp3 = n_Mult(tmp1, a, basecoeffs());
    tmp4 = n_Mult(tmp2, b, basecoeffs());
    n_InpAdd(tmp3, tmp4, basecoeffs());
    n_Delete(&tmp4, basecoeffs());

    n_InpMult(tmp1, c, basecoeffs());
    n_InpMult(tmp2, d, basecoeffs());
    n_InpAdd(tmp1, tmp2, basecoeffs());
    n_Delete(&tmp2, basecoeffs());

    set(i, j, tmp3);
    set(i, k, tmp1);
    n_Delete(&tmp1, basecoeffs());
    n_Delete(&tmp3, basecoeffs());
  }
}



//reduce all entries mod p. Does NOT change the coeffs type
void bigintmat::mod(number p, coeffs c)
{
  // produce the matrix in Z/pZ
  // CF: TODO rewrite using QuotRem and not the map
  coeffs coe = numbercoeffs(p, c);
  nMapFunc f1 = n_SetMap(basecoeffs(), coe);
  nMapFunc f2 = n_SetMap(coe, basecoeffs());
  number tmp1, tmp2;
  for (int i=1; i<=row; i++)
  {
    for (int j=1; j<=col; j++)
    {
      tmp1 = get(i, j);
      tmp2 = f1(tmp1, basecoeffs(), coe);
      n_Delete(&tmp1, basecoeffs());
      tmp1 = f2(tmp2, coe, basecoeffs());
      set(i, j, tmp1);
      n_Delete(&tmp1, basecoeffs());
      n_Delete(&tmp2, coe);
    }
  }
}

//coerce matrix to Z/pZ of type n_Zn
//NOT inplace! returns a copy
bigintmat* bigintmat::inpmod(number p, coeffs c)
{
  coeffs coe = numbercoeffs(p, c);
  return bimChangeCoeff(this, coe);
}

void bimMult(bigintmat *a, bigintmat *b, bigintmat *c)
{
  if (!nCoeffs_are_equal(a->basecoeffs(), b->basecoeffs())) {
    Werror("Error in bimMult. Coeffs do not agree!");
    return;
  }
  if ((a->rows() != c->rows()) || (b->cols() != c->cols()) || (a->cols() != b->rows())) {
    Werror("Error in bimMult. Dimensions do not agree!");
    return;
  }
  bigintmat *tmp = bimMult(a, b);
  c->copy(tmp);
  
  delete tmp;
}

//CF: TODO make work for Z/nZ, no inv.
number solvexA(bigintmat *A, bigintmat *b, bigintmat *x) {
  // Bestimme Pseudoinverse, mulipliziere von rechts an b und speichere Ergebnis in x.
  // Gebe Nenner von Pseudoinversen zurück.
  if (!nCoeffs_are_equal(A->basecoeffs(), b->basecoeffs())) {
    Werror("Error in solvexA. Coeffs do not agree!");
    return 0;
  }
  if (!nCoeffs_are_equal(A->basecoeffs(), x->basecoeffs())) {
    Werror("Error in solvexA. Coeffs do not agree!");
    return 0;
  }
  if ((x->rows() != b->rows()) || (x->cols() != A->rows()) || (A->cols() != b->cols())) {
    Werror("Error in solvexA. Dimensions do not agree!");
    return 0;
  }
  bigintmat *m = new bigintmat(A->rows(), A->cols(), A->basecoeffs());
  number t = A->pseudoinv(m);
  if (!(n_IsZero(t, A->basecoeffs()))) {
    bimMult(b,m,x);
    delete m;
    return t;
  }
  delete m;
  return 0;
}

//CF TODO make work for Z/nZ, don't use gauss
int kernbase (bigintmat *a, bigintmat *c, number p, coeffs q) {
  number temp;
  number minusone = n_Init(-1, a->basecoeffs());
  bigintmat *tmp = new bigintmat(a->rows(), 1, a->basecoeffs());
  bigintmat *m = a->modgauss(p, q);
  int dim = 0;
  for (int i=1; (i<=m->rows()) && (i <= m->cols()); i++) {
    temp = m->get(i, i);
    if (n_IsZero(temp, m->basecoeffs())) {
      dim++;
      m->set(i, i, minusone);
      m->getcol(i, tmp);
      c->setcol(dim, tmp);
    }
    n_Delete(&temp, m->basecoeffs());
  }
  n_Delete(&minusone, m->basecoeffs());
  delete tmp;  
  delete m;
  return dim;
}


//create Z/nA of type n_Zn
coeffs numbercoeffs(number n, coeffs c) {
  mpz_t p;
  number2mpz(n, c, p);
  ZnmInfo *pp = new ZnmInfo;
  pp->base = p;
  pp->exp = 1;
  coeffs nc = nInitChar(n_Zn, (void*)pp);
  return nc;
}

bool nCoeffs_are_equal(coeffs r, coeffs s) {
  if ((r == NULL) || (s == NULL))
    return false;
  if ((getCoeffType(r)==n_Z) && (getCoeffType(s)==n_Z))
    return true;
  if ((getCoeffType(r)==n_Zp) && (getCoeffType(s)==n_Zp)) {
    if (r->ch == s->ch)
      return true;
    else
      return false;
  }  
  // n_Zn stimmt wahrscheinlich noch nicht
  if ((getCoeffType(r)==n_Zn) && (getCoeffType(s)==n_Zn)) {
    if (r->ch == s->ch)
      return true;
    else
      return false;
  }
  if ((getCoeffType(r)==n_Q) && (getCoeffType(s)==n_Q))
    return true;
  // FALL n_Zn FEHLT NOCH!
    //if ((getCoeffType(r)==n_Zn) && (getCoeffType(s)==n_Zn))
  return false;
}

number bigintmat::content()
{
  coeffs r = basecoeffs();
  number g = get(1,1), h;
  int n=rows()*cols();
  for(int i=1; i<n && !n_IsOne(g, r); i++) {
    h = n_Gcd(g, view(i), r);
    n_Delete(&g, r);
    g=h;
  }
  return g;
}
void bigintmat::simplifyContentDen(number *d)
{
  coeffs r = basecoeffs();
  number g = n_Copy(*d, r), h;
  int n=rows()*cols();
  for(int i=0; i<n && !n_IsOne(g, r); i++) {
    h = n_Gcd(g, view(i), r);
    n_Delete(&g, r);
    g=h;
  }
  *d = n_Div(*d, g, r);
  if (!n_IsOne(g, r))
  skaldiv(g);
}
