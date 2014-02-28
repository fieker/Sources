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
void nforder::init() {
  // Gibt es eine Multtable, so gibt es keine Baseorder
  baseorder = NULL;
  basis = NULL;
  // Discriminante wird erst berechnet, wenn sie benötigt wird
  discriminant = NULL;
  divisor = NULL;
  flags = 0;
  multtable = NULL;
  m_coeffs = NULL;
  setOneIsOne();
}

nforder::nforder(int dim, bigintmat **m,const coeffs q) {
  init();
  m_coeffs = q;
  dimension = dim;
  multtable = (bigintmat**)(omAlloc(dim*sizeof(bigintmat*)));
  for (int i=0; i<dim; i++) {
    multtable[i] = new bigintmat(m[i]);
  }
}

nforder::nforder(nforder *o, bigintmat *base, number div, const coeffs q) {
   init();
   m_coeffs = q;  
   basis = new bigintmat(base);
   //neue Ordnung erzeugen und übergebene Daten kopieren
   baseorder = new nforder(o, 1); 
   //Gibt es eine Baseorder, brauchen wir keine Multtable. Könnte aber evtl. generiert werden
   multtable = NULL;
   divisor = n_Copy(div,basecoeffs());
   dimension = o->getDim();
   discriminant = NULL;
}

nforder::nforder(nforder *o, int) {
  init();
  m_coeffs = o->basecoeffs();
  ::Print("copy called: %lx\n", m_coeffs);
  // Kopiert die Daten der übergebenen Ordnung auf die erzeugte
  if (o->discriminant)
    discriminant = n_Copy(o->discriminant,  basecoeffs());
  // get-Funktionen liefern immer nur Kopien der Attribute zurück
  dimension = o->getDim();
  multtable = (bigintmat **)omAlloc(dimension*sizeof(bigintmat*));
  if (!o->getMult(multtable)) {
    omFree(multtable);
    multtable = NULL;
  }
  baseorder = o->getBase();
  basis = o->getBasis();
  if (o->divisor)
    divisor = n_Copy(o->divisor,  basecoeffs());
}

void nforder::Print() {
  ::Print("Order:\nof dimension %d\n", dimension);
  if (discriminant && !n_IsZero(discriminant, m_coeffs)) {
    ::Print("and discriminant: ");
    StringSetS("");
    n_Write(discriminant, m_coeffs);
    ::Print("%s\n", StringEndS());
  }
//  coeffs
  if (multtable) {
    ::Print("Multiplication table:\n");
    for(int i=0; i<dimension; i++) {
      char * m = multtable[i]->String();
      ::Print("%d: %s\n", i, m);
    }
  }

  if (baseorder) {
    ::Print("as extension of: ");
    baseorder->Print();
    ::Print("with basis:\n");
    char * m = basis->String();
    ::Print("%s", m);
    ::Print("and denominator: ");
    StringSetS("");
    n_Write(divisor, m_coeffs);
    ::Print("%s\n", StringEndS());
  }

  ::Print("Flags: %lx\n", flags);
}

// Was bei Desktruktor ändern?
nforder::~nforder() {
  if (multtable != NULL) {
    // Falls es eine Multtable gab, werden zunächst sämtliche n Matrizen gelöscht, und anschließend das (durch malloc) erzeugte Array of Matrix
    for (int i=0; i<dimension; i++)
      delete multtable[i];
    omFree(multtable);
  }
  else
  {
    // andernfalls werden baseorder und basis gelöscht
    delete baseorder;
    delete basis;
    if (divisor) n_Delete(&divisor, basecoeffs());
  }
  if (discriminant) n_Delete(&discriminant, basecoeffs());
}

/*
////////////////////////////////////
//////// Private Funktionen ////////
///////////////////////////////////*/
/*_______________-1_______________ */
void nforder::calcdisc() {
  // Determinante von Spurmatrix ist die Discriminante
  puts("calcdisc called");
  if (baseorder == NULL) {
    bigintmat *m = traceMatrix();
    discriminant = m->det();
    delete m;
  }
  else
  {
    puts("not easy");
    number prod = n_Init(1, basecoeffs());
    number tmp, tmp2;
    for (int i=1; i<=dimension; i++) {
      tmp2 = basis->get(i, i);
      tmp = n_Mult(prod, tmp2, basecoeffs());
      n_Delete(&tmp2, basis->basecoeffs());
      n_Delete(&prod, basecoeffs());
      prod = tmp;
    }
    number disc = baseorder->getDisc();
    number detquad = n_Mult(prod, prod, basis->basecoeffs());
    discriminant = n_Mult(disc, detquad, basecoeffs());

    for (int i=1; i<=2*dimension; i++) {
      tmp = n_Div(discriminant, divisor, basecoeffs());
      n_Delete(&discriminant, basecoeffs());
      discriminant = tmp;
    }
    n_Delete(&detquad, basis->basecoeffs());
    n_Delete(&disc, baseorder->basecoeffs());
  }
}

bigintmat *nforder::traceMatrix() {
  bigintmat *m = new bigintmat(dimension, dimension, basecoeffs());
  bigintmat *base1 = new bigintmat(1, dimension, basecoeffs());
  bigintmat *base2 = new bigintmat(1, dimension, basecoeffs());
  bigintmat *mm = new bigintmat(dimension, dimension, basecoeffs());
  number sum;
  number temp;
  number t1;
  
  for (int i=1; i<=dimension; i++) {
    for (int j=1; j<=dimension; j++) {
      // Berechnet Produkt von Basiselementen i und j und speichert es in base1
      makebase(base1, i);
      makebase(base2, j);
      elMult(base1, base2);
      // Schreibt Abbildungsmatrix der Multiplikation mit base1 in mm
      multmap(base1, mm);
      sum = n_Init(0,basecoeffs());
      // Berechnet Spur der Abbildungsmatrix...
      for (int k=1; k<=dimension; k++) {
        temp = n_Copy(sum, basecoeffs());
        t1 = mm->get(k,k);
        n_Delete(&sum, basecoeffs());
        sum = n_Add(t1, temp, basecoeffs());
        n_Delete(&t1, basecoeffs());
        n_Delete(&temp, basecoeffs());
      }
      // ... und schreibt diese als Eintrag in Matrix m an Stelle (i, j)
      m->set(i, j, sum, basecoeffs());
      n_Delete(&sum, basecoeffs());
    }
  }
  delete base1;
  delete base2;
  delete mm;
                          // Hier stand vorher: n_Delete(&t1, basecoeffs()); Funktion noch testen.Funktion gestestet. Müsste funktionieren.
  return m;
}
////////////////////////////////////
////// Öffentliche Funktionen //////
////////////////////////////////////

/*_____________+1_______________ */
number nforder::getDisc() {
  // Falls Discriminante bisher noch nicht berechnet wurde, berechne diese
  if (!discriminant || n_IsZero(discriminant, basecoeffs())) {
    calcdisc();
  }
  return n_Copy(discriminant, basecoeffs());
}

int nforder::getDim() {
  return dimension;
}

bigintmat *nforder::getBasis() {
  // Falls basis ein NULL-Pointer ist, liefere NULL zurück, andernfalls liefere eine Kopie von basis
  if (basis == NULL)
    return NULL;
  bigintmat *m = new bigintmat(basis); //wenn Fehler dann hier
  return m;
}

bool nforder::getMult(bigintmat **m) {
  // Falls multtable ein NULL-Pointer ist, liefere NULL zurück, andernfalls erzeuge neues Array of Matrix, kopiere die Matrizen aus multtable dort hinein, und gib es zurück
  if (multtable == NULL) {
    return false;
  }
  for (int i=0; i<dimension; i++)
  {
    m[i] = new bigintmat(multtable[i]);
  }
  return true;
}


number nforder::getDiv() {
  return n_Copy(divisor, basecoeffs());
}

nforder *nforder::getBase() {
  // Falls baseorder ein NULL-Pointer ist, liefere NULL zurück, andernfalls liefere eine Kopie von baseorder
  if (baseorder == NULL)
    return NULL;
  nforder *o = new nforder(baseorder, TRUE);
  return o;
}

number nforder::foundBasis(bigintmat *b) {
  // Falls es keine baseorder gibt, liefere als Nenner 0 zurück
  if (baseorder == NULL) {
    return n_Init(0, basecoeffs());
  }
  bigintmat *check = baseorder->getBasis();
  if (check == NULL) {
    // Falls baseorder keine Basis hat (und damit also eine multtable), sind die eigene Basis und Nenner schon die gewünschten
    delete b;
    b=getBasis();
    return divisor;
  }
  
  // Andernfalls greife auf die foundBasis-Funktion von baseorder zu (speichere Nenner in di, und Basis in fb).
  bigintmat *fb = new bigintmat(dimension, dimension,basecoeffs());
  number di = baseorder->foundBasis(fb);
  bigintmat *tmp = new bigintmat(1, dimension,basecoeffs());
  bigintmat *sum = new bigintmat(1, dimension,basecoeffs());
  for (int i=1; i<=dimension; i++) {
    // Laufe durch alle Basiselemente durch
    sum->zero();
    for (int j=1; j<=dimension; j++) {
      // Laufe durch die Basiselemente von baseorder und speichere die Darstellung des j-ten Basiselements in tmp
      fb->getrow(j, tmp);
      // Multipliziere sie mit basis[i, j] (Koeff. von baseorder-Basiselem. j bei Linearkombination von this-Basiselement i durch Basiselemente von baseorder)
      tmp->skalmult(basis->get(i, j), basis->basecoeffs());
      // Addiere diese Darstellung für alle j auf
      sum->add(tmp);
    }
    // Das Ergebnis ist die grundlegende Basisdarstellung (Lin.-Komb. von Basiselementen der untersten Ordnung) von this-Basiselement i
    b->setrow(i, sum);
  }
  
  delete check;
  delete tmp;
  delete sum;
  delete fb;
  // Der Nenner ist der grundlegende Nenner von baseorder multipliziert mit this.Nenner
  number t = n_Mult(divisor,di,basecoeffs());
  n_Delete(&di,basecoeffs());
  return t;
}

// _____________+2_______________ 
void nforder::elAdd(bigintmat *a, bigintmat *b) {
  if ((a->rows() != 1) || (a->cols() != dimension) || (b->rows() != 1) || (b->cols() != dimension)) {
    // Kein Zeilenvektor der korrekten Größe
    Werror("Error in elAdd");
  }
  else {
    a->add(b);
  }
}

void nforder::elSub(bigintmat *a, bigintmat *b) {
  if ((a->rows() != 1) || (a->cols() != dimension) || (b->rows() != 1) || (b->cols() != dimension)) {
    // Kein Zeilenvektor der korrekten Größe
    Werror("Error in elSub");
  }
  else {
    a->sub(b);
  }
}

void nforder::elMult(bigintmat *a, bigintmat *b) {
  if ((a->rows() != 1) || (a->cols() != dimension) || (b->rows() != 1) || (b->cols() != dimension)) {
    // Kein Zeilenvektor der korrekten Größe
    Werror("Error in elMult");
  }
  else {
    // Teil mit multtable funktioniert
    if (multtable != NULL) {
      // Multiplikation mit Hilfe von Multiplikationstabelle
      // Zu Grunde liegende Formel: Basis w_i; Für alpha = sum a_i*w_i und beta = sum b_i*w_i gilt:
      // alpha*beta = sum sum a_i*b_j*w_i*w_j
      bigintmat *sum = new bigintmat(1, dimension, a->basecoeffs());
      bigintmat *tmp = new bigintmat(1, dimension, a->basecoeffs());
      number ntmp;
      for (int i=1; i<=dimension; i++) {
        // Laufe mit i durch Basiselemente
        for (int j=1; j<=dimension; j++) {
          // Laufe mit j durch Basiselemente
          // Speichere Produkt von Basiselem. i mit Basiselem. j als Koeff.vektor in tmp
          
          multtable[i-1]->getrow(j, tmp);
          // Multipliziere ihn mit a[i] und b[j]
          ntmp = n_Mult(a->get(i-1), b->get(j-1), a->basecoeffs());
          tmp->skalmult(ntmp, a->basecoeffs());
          
          n_Delete(&ntmp, a->basecoeffs());
          // und addiere alles auf
          sum->add(tmp);
        }
      }
      delete tmp;
      // Am Ende überschreibe a mit dem Ergebnis
      for (int i=0; i<dimension; i++)
        a->set(i, sum->get(i));
      delete sum;
    }
    else {
    // Multiplikation mit hilfe von baseorder:
      bigintmat *suma, *sumb;
      suma = new bigintmat(1, dimension, a->basecoeffs());
      sumb = new bigintmat(1, dimension, a->basecoeffs());
    // Produkt von a (b) mit basis liefert Koeff-Vektor von a*divisor (b*divisor) in baseorder
      bimMult(a, basis, suma);
      bimMult(b, basis, sumb);
      // Multipliziere Elemente in baseorder (und speichere in suma)
      baseorder->elMult(suma, sumb);
      // Dividiere einmal den Nenner heraus, um Koeff-Vektor von a*b*divisor in baseorder zu erhalten
      suma->skaldiv(divisor);
      // Löse das Gleichungssysten x*A = b, wobei A = basis, und b = Koeffvektor von a*b*divisor in baseorder.
      // Die Lösung wird durch den Nenner der Lösung (Rückgabewert von solvexA) geteilt (geht hier immer auf, da wir eine Basis haben)
      // und enthält dann den Koeff-Vektor von a*b in this (Ergebnis wird wieder in a gespeichert)
      a->skaldiv(solvexA(basis, suma, a));
      delete suma;
      delete sumb;
    }
  }
}


void nforder::multmap(bigintmat *a, bigintmat *m) {
  if ((m->cols() != dimension) || (m->rows() != dimension)) {
    Werror("Error in multmap");
    return;
  }
  bigintmat *bas = new bigintmat(1, dimension, basecoeffs());
  for (int i=1; i<=dimension; i++) {
    // Durchläuft alle Basiselemente
    // Multipliziert i-tes Basiselement mit a
    makebase(bas, i);
    elMult(bas, a);
    // Schreibt Ergebnis in i-te Zeile der Matrix m. Am Ende ist m dann die Abbildungsmatrix der Multiplikation mit a
    m->setrow(i, bas);
  }
  delete bas;
}

/*________________1_______________ */
void makebase(bigintmat *m, int i) {
  if (((m->rows() == 1) && (i <= m->cols())) || ((m->cols() == 1) && (i <= m->rows()))) {
    // Falls m Zeilen- oder Spaltenvektor ist, setze alle Einträge auf 0 und Eintrag i auf 1 (Koeff-Vektor des i-ten Basiselements)
    number t1 = n_Init(0,m->basecoeffs());
    for (int j=0; ((j<m->rows()) || (j<m->cols())); j++) {
      m->set(j, t1);
      
    }
    n_Delete(&t1,m->basecoeffs());
    number t2 = n_Init(1,m->basecoeffs());
    m->set(i-1, t2);
    n_Delete(&t2,m->basecoeffs());
  }
  else
    Werror("Error in makebase. Not a vector.");
}

////////////////////////////////////
//////////// 2 Round 2 /////////////
////////////////////////////////////

bigintmat *radicalmodpbase(nforder *o, number p, coeffs c) {
  
  number dimen = n_Init(o->getDim(), o->basecoeffs());
  int n = o->getDim();

  bigintmat *m, *bas;
  // Berechnet F_p-Basis von I_p/pI_p (Radical mod p)
  // Dazu: 
  if (n_Greater(p, dimen, c)) {
    // Falls Primzahl größer gleich Dimension der Ordnung, so berechne Kern der Spurmatrix modulo p.
    m = o->traceMatrix();
    bas = new bigintmat(n, 1, o->basecoeffs());
  }
  else {
    // Sonst: Berechne Kern der Abbildung x -> x^(p^j) mod p, wobei j>0 mit p^j >= dimension
    int j = 1;
    // ex als number, oder reicht long long int?
    // Finde j von oben und berechne p^j
    number ex = n_Init(1, o->basecoeffs());
    number temp;
    while (n_Greater(dimen, ex, o->basecoeffs())) {
      temp = n_Mult(ex, p, o->basecoeffs());
      n_Delete(&ex, o->basecoeffs());
      ex = temp;
      j++;
    }
    
    // Berechne Abbildungsmatrix der oben genannten Abbildung und speichere diese in m (genauere Erklärung dazu: Siehe multmap())
    m = new bigintmat(n, n, o->basecoeffs());
    bas = new bigintmat(1, n, o->basecoeffs());
    bigintmat *prod = new bigintmat(1, n, o->basecoeffs());
    
    number klauf;
    number eins = n_Init(1, o->basecoeffs());
    
    for (int i=1; i<=n; i++) {
      makebase(bas, i);
      prod->copy(bas);
      klauf = n_Init(1, o->basecoeffs());
      for (; n_Greater(ex, klauf, o->basecoeffs());) {
        o->elMult(prod, bas);
        prod->mod(p, c);
        temp = n_Add(klauf, eins, o->basecoeffs());
        n_Delete(&klauf, o->basecoeffs());
        klauf = temp;
      }
      n_Delete(&klauf, o->basecoeffs());
      m->setcol(i, prod);
    }

    delete prod;
    n_Delete(&ex, o->basecoeffs());
    n_Delete(&eins, o->basecoeffs());
    
  }
  
  bigintmat *kbase = new bigintmat(n, n, o->basecoeffs());

  
  // Speichere Basiselemente von Kern der Matrix m (Spurmatrix oder Abbildungsmatrix, je nach if-else-Fall) (von Z/pZ -> Z/pZ) in kbase (ersten kdim Spalten bilden Basis)
  int kdim = kernbase(m, kbase, p, c);
  // Schreibe für jedes i=1,, .., dimension p*(i-tes Basiselement) als Spalten in Matrix gen, dahinter die oben errechnete Basis vom Kern
  // Wir erhalten (als Spalten) ein Erzeugendensystem vom Kern von Z->Z/pZ: x->x^(p^j)
  bigintmat *gen = new bigintmat(n, n+kdim, o->basecoeffs());
  
  for (int i=1; i<=n; i++) {
    makebase(bas, i);
    bas->skalmult(p, c);
    gen->setcol(i, bas);
  }
  for (int i=1; i<=kdim; i++) {
    kbase->getcol(i, bas);
    gen->setcol(i+n, bas);
  }
  
  // HNF auf EZS anwenden liefert (als letzten dimension Spalten) eine Basis des Kerns
  gen->hnf();
  bigintmat *tmp = new bigintmat(n, 1, o->basecoeffs());
  bigintmat *nbase = new bigintmat(n, n, o->basecoeffs());
  // Schreibe diese als Spalten in nbase und gib nbase zurück
  for (int i=1; i<=n; i++) {
    gen->getcol(gen->cols()-n+i, tmp);
    nbase->setcol(i, tmp);
  }
  
  n_Delete(&dimen, o->basecoeffs());
  delete m;
  delete bas; 
  delete kbase;
  delete gen;
  delete tmp;
  return nbase;

}

number multring(bigintmat *nbase, nforder *o, number p) {
  number divi;
  int n = o->getDim();
  bigintmat *inv = new bigintmat(n, n, o->basecoeffs());
  divi = nbase->pseudoinv(inv);
  // Zusammenbau der "langen" Matrix
  bigintmat *lon = new bigintmat(n, n, o->basecoeffs());
  bigintmat *oldlon;
  bigintmat *mm = new bigintmat(n, n, o->basecoeffs());
  bigintmat *temp = new bigintmat(1, n, o->basecoeffs());
  bigintmat *nochnetemp = new bigintmat(n, n, o->basecoeffs());

  for (int i=1; i<=n; i++) {
    nbase->getrow(i, temp);
    o->multmap(temp, mm);
    bimMult(mm, inv, nochnetemp);
    mm->copy(nochnetemp);
    oldlon = new bigintmat(lon);
    delete lon;
    lon = new bigintmat(n, (i+1)*n, o->basecoeffs());
    lon->concatcol(oldlon, mm);
    delete oldlon;
  }
  
  lon->skaldiv(divi);

  number p2, zwei;
  zwei = n_Init(2, o->basecoeffs());
  p2 = n_Mult(p, zwei, o->basecoeffs());
  lon->modhnf2(p2, o->basecoeffs());
  
  n_Delete(&zwei, o->basecoeffs());
  n_Delete(&p2, o->basecoeffs());
  bigintmat *red = new bigintmat(n, n, o->basecoeffs());
  bigintmat *ttemp = new bigintmat(n, 1, o->basecoeffs());
  
  for (int i=1; i<=n; i++) {
    lon->getcol(n*n+i, ttemp);
    red->setcol(i, ttemp);
  }
  number divisor = red->pseudoinv(nbase);
  
  delete inv;
  delete lon;
  delete mm;
  delete temp;
  delete red;
  delete ttemp;
  delete nochnetemp;
  n_Delete(&divi, o->basecoeffs());
  return divisor;
}

nforder *onestep(nforder *o, number p, coeffs c) {
  // Berechne F_p-Basis von I_p/pI_p
  bigintmat *basis;
  basis = radicalmodpbase(o, p, c);
  // Wir benötigen die Matrix Zeilen- statt Spaltenweise
  bigintmat *tmp = basis->transpose();
  delete basis;
  basis = tmp;
  // Bestimme Basis vom Ring der Multiplikatoren (speicher diese in basis), und Nenner davon (in divisor)
  number divisor = multring(basis, o, p);
  // Erzeuge neue Ordnung, der o zu Grunde liegt, mit Basis basis und Nenner divisor
  nforder *no = new nforder(o, basis, divisor, o->basecoeffs());
  
  
  delete basis;
  n_Delete(&divisor, o->basecoeffs());
  return no;
}

nforder *pmaximal(nforder *o, number p) {
  coeffs c = o->basecoeffs();
  number disc = o->getDisc();
  number olddisc = n_Init(0, c);
  nforder *no;
  nforder *otemp = new nforder(o, 1);
  // Solange onestep etwas tut (also die Diskriminante sich verändert) führe erneut onestep aus
  while (!n_Equal(olddisc, disc, otemp->basecoeffs())) {
    n_Delete(&olddisc, otemp->basecoeffs());
    olddisc = disc;
    no = onestep(otemp, p, c);
    delete otemp;
    otemp = no;
    disc = otemp->getDisc();
  }
  n_Delete(&disc, otemp->basecoeffs());
  n_Delete(&olddisc, otemp->basecoeffs());
  return no;
}

/*
// Zum Round2 fehlt noch die Faktorisierung der Diskriminante. Daher auch noch nicht getestet
nforder *round2(nforder *o) {
  nforder *otemp = new nforder(o,basecoeffs());
  number p = otemp->getsmallestsqprime(); // Benötigt kleinste Primzahl, die die Disc. quadratisch teilt
  nforder *no;
  number eins = n_Init(1, basecoeffs());
  number tmp;
  while (n_GreaterZero(p,basecoeffs())) {
    // Laufe durch Primzahlen p, die die Ordnung quadratisch teilen, und erzeuge p-maximale Ordnung
    no = pmaximal(otemp, p);
    delete otemp;
    otemp = no;
    // Nimm nächstgrößere Primzahl, welche die Ordnung quadratisch teilt
    tmp = n_Add(p,eins, basecoeffs());
    p = otemp->getsmallestsqprime(tmp); // Benötigt kleinste Primzahl größer tmp, die die Disc. quad. teilt
    n_Delete(&tmp, basecoeffs());
  }
  n_Delete(&p, basecoeffs());
  n_Delete(&eins, basecoeffs());
  return otemp;
}
*/


void nforder::createmulttable(bigintmat **a) {
  // Falls es eine Multtable gibt, liefere eine Kopie davon zurück
  if (multtable != NULL) {
    for (int i=0; i<dimension; i++) {
      a[i] = new bigintmat(multtable[i]);
    }
  }
  else {
    // Sonst berechne sie auf kanonische Art und Weise
    bigintmat *bas = new bigintmat(1, dimension, basecoeffs());
    for (int i=0; i<dimension; i++) {
      makebase(bas, i+1);
      a[i] = new bigintmat(dimension, dimension, basecoeffs());
      multmap(bas, a[i]);
    }
  }
}
