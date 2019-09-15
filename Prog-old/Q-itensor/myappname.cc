#include "itensor/all.h"
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <typeinfo>
using namespace itensor;
using namespace std;

int main()
{
auto q = QN("Nb=",2);
auto q1 = QN("Nb=",3);
auto q2 = QN("Nb=",4);

//auto b = QN({0,3});
//auto b1 = QN({1,3});
//auto b2 = QN({2,3});
auto g = QN({2,1});

auto b = QN({0,3});
auto b1 = QN({2,3});
//auto b2 = QN({2,3});



printFull(b+b1);
//printFull(b1+b2);
printFull(g+g);
//Arrow dir = Out;
//auto I = IQIndex("I",Index("b1",2),b,Index("b2",2),b1,Index("b3",2),b2,Out);
//auto J = IQIndex("I",Index("j1",2),b,Index("j2",2),b1,Index("j3",2),b2,Out);

auto I = IQIndex("I",Index("b1",1),b,Index("b2",1),b1,Out);
auto J = IQIndex("J",Index("j1",1),b,Index("j2",1),b1,Out);


Print(I);
Print(I.m()); 
Print(I.nblock()); 
Print(I.index(2)); 
Print(I.qn(2));    
Print(I.dir()); 
Print(dag(I)); 

auto A = IQTensor(I,dag(I));
auto B = randomTensor(QN({0,3}),I,J,prime(dag(I)), prime(dag(J)) );
//A.randomTensor();
PrintData(A);
Print(B);
PrintData(B);

Print(div(B));

IQTensor U(I,J),S,V;
svd(B,U,S,V,{"Cutoff",1E-4});
Print(U);

Print(sqr(norm(B-U*S*V)/norm(B))); //prints: 1E-4

//for(auto& iq : I)
//  {
//  println(iq.index);
//  println(iq.qn);
//  }




auto v = stdx::reserve_vector<IndexQN>(5);
v.emplace_back(Index("I+2",4),QN(+2));
v.emplace_back(Index("I+1",8),QN(+1));
v.emplace_back(Index("I_0",10),QN(0));
v.emplace_back(Index("I-1",8),QN(-1));
v.emplace_back(Index("I-2",4),QN(-2));

auto I1 = IQIndex("I",std::move(v),Out,0);



//Print(b.mod(1)); //prints: b.mod(1) = 3
//Print(b.mod(2)); //prints: b.mod(2) = 2


//println(Nb(q));
//printFull(q);

//println(q);
//Print(q[0]);
//Print(q[1]);

//Print(q.mod(1));

return 0;
}



