//
// test script for ord command
//
pagelength = 10000;
ring r1 = 32003,(x,y,z),wp(2,3,2);
r1;
poly s1=x2+y2+z2;
s1;
ord(s1);
"--------------------------------";
ring r2=23,(x,y,z),ws(2,3,2);
r2;
poly s1=x2+y2+z2;
s1;
ord(s1);
"-------------------------------";
ring r3=0,(x,y,z(1..18)),ws(4,2,3,4,5,2,3,1,5,3,2,2,2,3,3,3,4,5,3,1);
r3;
poly s2=x2+y2+z(1)^8+z(12)^18;
s2;
ord(s2);
"------------------------------";
vector v=[x3,y2,z(15)^2,z(2)^2];
v;
ord(v);
"--------------------------------";
listvar(all);
kill r1,r2,r3;
LIB "tst.lib";tst_status(1);$;
