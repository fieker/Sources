//////////////////////////////////////////
//////////////////////////////////////////
////     ideals in oforder    ////////////
//////////////////////////////////////////
//////////////////////////////////////////
#ifndef NFORDER_IDEAL_HPP
#define NFORDER_IDEAL_HPP

class nforder_ideal
{
private:
  ////////////////////////////////////
  ////////// Membervariablen /////////
  ////////////////////////////////////
  number norm, norm_den, min, min_den;
  coeffs ord;  // but of dynamic type order! (via the black-box-stuff)
  bigintmat *basis; 
  
  
public:
  
  ////////////////////////////////////
  /// 0 Konstruktoren/Destruktoren ///
  ////////////////////////////////////
   //Lädt entweder Multiplikationstabelle von location, oder legt Ordnung o zu Grunde (mit Basis base und Hauptnenner div)  
  nforder_ideal();
  nforder_ideal(bigintmat *basis, 
                   number * min, number * min_den, 
                   number * norm, number * norm_den, 
                   const coeffs ord); 
  nforder_ideal(bigintmat *basis, const coeffs * q); 
  nforder(nforder_ideal *i);
  
  ~nforder_ideal();
  void Print();

  ////////////////////////////////////
  // +1 Zugriff auf Membervariablen //
  ////////////////////////////////////
  
  number getNorm(), getNormDen(), getMin(), getMinDen();
  inline coeffs order() const { return ord; }
  bigintmat *getBasis();
  
  
  ////////////////////////////////////
  ////// +2 Elementoperationen ///////
  ////////////////////////////////////
  // Addiert/Subtrahiert/Multipliziert zu a das Element b hinzu  
  void elAdd(nforder_ideal *a, nforder_ideal *b);
  void elMult(nforder_ideal *a, nforder_ideal *b);
  void elDiv(nforder_ideal *a, nforder_ideal *b);
  void elMeet(nforder_ideal *a, nforder_ideal *b);
};

#endif
