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
  coeffs ord;  // but of dynamic type order! (as cring)
  bigintmat *basis; 
  
  
public:
  
  ////////////////////////////////////
  /// 0 Konstruktoren/Destruktoren ///
  ////////////////////////////////////
  nforder_ideal();
  void init();
  nforder_ideal(bigintmat *basis, 
                   number * min, number * min_den, 
                   number * norm, number * norm_den, 
                   const coeffs ord); 
  nforder_ideal(bigintmat *basis, const coeffs q); 
  nforder_ideal(nforder_ideal *i, int);
  
  ~nforder_ideal();
  void Write();
  char * String();
  void Print();

  ////////////////////////////////////
  // +1 Zugriff auf Membervariablen //
  ////////////////////////////////////
  
  number getNorm(), getNormDen(), getMin(), getMinDen();
  inline coeffs order() const { return ord; }
  inline bigintmat * viewBasis() {return basis;};
  inline void setMinTransfer(number a, number b){min = a; min_den = b;}
  inline void setNormTransfer(number a, number b){norm = a; norm_den = b;}
  
  ////////////////////////////////////
  ////// +2 Elementoperationen ///////
  ////////////////////////////////////
  // Addiert/Subtrahiert/Multipliziert zu a das Element b hinzu  
};

nforder_ideal* nf_idAdd(nforder_ideal *a, nforder_ideal *b);
nforder_ideal* nf_idMult(nforder_ideal *a, nforder_ideal *b);
nforder_ideal* nf_idDiv(nforder_ideal *a, nforder_ideal *b);
nforder_ideal* nf_idMeet(nforder_ideal *a, nforder_ideal *b);
#endif
