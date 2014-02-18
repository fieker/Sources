//////////////////////////////////////////
//////////////////////////////////////////
////     Einfache Ordnungs-Klasse     ////
//////////////////////////////////////////
// Kira Kraft, Jan Albert, Marco Sieben //
/////// Praktikum bei Herrn Fieker ///////
//////////////////////////////////////////
//////////////////////////////////////////
/*
 * - Einer Ordnung liegt entweder eine Multiplikationstabelle zu Grunde, oder sie wird durch eine O-Basis erzeugt, wobei O eine andere Ordnung ist
 * - Ordnungselemente werden als Koeffizientenvektoren (Klasse matrix, Zeilenvektoren!) dargestellt und können addiert/subtrahiert/multipliziert werden
 * - Die Diskriminante kann von selbst berechnet werden
*/
#ifndef NFORDER_HPP
#define NFORDER_HPP

class nforder
{
private:
  ////////////////////////////////////
  ////////// Membervariablen /////////
  ////////////////////////////////////
  number discriminant;
  int dimension;
  coeffs m_coeffs;
  bigintmat **multtable; // Multiplikationstabelle als Array von Matrizen ... 
  nforder *baseorder; // ... oder zu Grunde liegende Ordnung  
  bigintmat *basis; // Lin.Komb. der Basiselemente von baseorder zu Basiselementen der Ordnung (Eine Zeile ^= ein Basiselement)  
  number divisor; // Hauptnenner der Linearkombination der Basiselemente  
    // Entweder multtable oder baseorder zeigt auf NULL - je nachdem, wie die Ordnung konstruiert wurde  
  
  ////////////////////////////////////
  /////// -1 Memberfunktionen ////////
  ////////////////////////////////////
  // Genauere Beschreibung aller Funktionen in der Funktionen.odt
  
public:
  void calcdisc(); // Berechnet Diskriminante
  
  
  ////////////////////////////////////
  /// 0 Konstruktoren/Destruktoren ///
  ////////////////////////////////////
   //Lädt entweder Multiplikationstabelle von location, oder legt Ordnung o zu Grunde (mit Basis base und Hauptnenner div)  
  nforder(int dim, bigintmat **m, const coeffs q); // (keine Übergabe von const char *, sondern von bigintmat *, diese wird kopiert und als multtable verwendet)
  nforder(nforder *o, bigintmat *base, number div, const coeffs q);
  nforder(nforder *o);
  
  ~nforder();
  void Print();

  ////////////////////////////////////
  // +1 Zugriff auf Membervariablen //
  ////////////////////////////////////
  
  number getDisc();
  int getDim();
  inline coeffs basecoeffs() const { return m_coeffs; }
  number getDiv();
  // Liefert Zeiger auf Kopier der Objekte zurück  
  bool getMult(bigintmat **m);
  nforder *getBase();
  bigintmat *getBasis();
  number foundBasis(bigintmat *b); // Liefert Basis in Abhängigkeit von zu unterst liegenden Ordnung
  
  
  
  ////////////////////////////////////
  ////// +2 Elementoperationen ///////
  ////////////////////////////////////
  // Addiert/Subtrahiert/Multipliziert zu a das Element b hinzu  
  void elAdd(bigintmat *a, bigintmat *b);
  void elSub(bigintmat *a, bigintmat *b);
  void elMult(bigintmat *a, bigintmat *b);
  
  ////////////////////////////////////
  //// +3 Funktionen für Round 2 /////
  ////////////////////////////////////
  // long long int getsmallestsqprime();
  /* Liefert kleinste Primzahl >= p, die die Diskriminante quadratisch teilt  */
  void multmap(bigintmat *a, bigintmat *m);
  bigintmat *traceMatrix();
  void createmulttable(bigintmat **a);
  
};

////////////////////////////////////
////// 1 Komfortfunktionen /////////
////////////////////////////////////
/* Setzt Vektor m auf (0,...,0,1,0,...,0) (i-ten Basisvektor)  */
void makebase(bigintmat *m, int i);

////////////////////////////////////
//////////// 2 Round 2 /////////////
////////////////////////////////////
/* Liefert bzgl. Primzahl p um eines größere Ordnung von o zurück */ 
nforder *onestep(nforder *o, number p, coeffs c);
/* Macht liefert p-maximale Ordnung von o zurück  */
nforder *pmaximal(nforder *o, number p, coeffs c);
/* Liefert Maximalordnung, ausgehend von o, zurück  */
nforder *round2(nforder *o); // Benötigt Faktorisierung der Diskriminanten
/* Liefert Basis von I_p(O)/pI_p(O)  */
bigintmat *radicalmodpbase(nforder *o, number p, coeffs c);
/* Berechnet die Basis mit Hilfe der langen Matrix  */
number multring(bigintmat* nbase, nforder *o, number p);

#endif
