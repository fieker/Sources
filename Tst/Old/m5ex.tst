  LIB "sing.lib";
  ring R=32003,(x,y,z),ds;
  // ---------------------------------------
  // hypersurface case (from series T[p,q,r]):
  int p,q,r = 3,3,4;
  poly f = x^p+y^q+z^r+xyz;
  tjurina(f);
  kbase(Tjurina(f));
  // Tjurina number = 8
  // ---------------------------------------
  // complete intersection case (from series P[k,l]):
  int k,l =3,2;
  ideal j=xy,x^k+y^l+z2;
  dim(std(j));          // Krull dimension
  size(minbase(j));     // minimal number of generators
  tjurina(j);           // Tjurina number
  module T=Tjurina(j);
  kbase(T);             // a sparse output of the k-basis of T1
  print(kbase(T));      // columns of matrix are a k-basis of T1
  // ---------------------------------------
  // general case (cone over rational normal curve of degree 4):
  ring r1=0,(x,y,z,u,v),ds;
  matrix m[2][4]=x,y,z,u,y,z,u,v;
  ideal i=minor(m,2);   // 2x2 minors of matrix m
  module M=T1(i);       // a presentation matrix of T1
  vdim(M);              // Tjurina number
  hilb(M);              // display of both Hilbert series
  intvec v1=hilb(M,1);  // first Hilbert series as intvec
  intvec v2=hilb(M,2);  // second Hilbert series as intvec
  v1;
  v2;
  v1[3];                // 3-rd coefficient of the 1-st Hilbert series
  module N=T2(i);
LIB "tst.lib";tst_status(1);$
