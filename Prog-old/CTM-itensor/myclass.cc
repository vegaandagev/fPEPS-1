#include "Header.h"
auto   Transpose_itensor( ITensor &, Index&, Index&  )  -> ITensor;


auto Init_Tensors(int& chi, int& D) -> vector<ITensor>
{
vector<ITensor> Ten_Iten;

auto Ain=Index("A_in", D, Atype, 0);
auto Bout=Index("B_out", D, Btype, 0);

Ten_Iten.push_back(randomTensor(Ain,prime(Ain,1),Bout,prime(Bout,1)));
Ten_Iten.push_back(randomTensor(Ain,prime(Ain,1),Bout,prime(Bout,1)));
Ten_Iten.push_back(randomTensor(Ain,prime(Ain,1),Bout,prime(Bout,1)));
Ten_Iten.push_back(randomTensor(Ain,prime(Ain,1),Bout,prime(Bout,1)));
Ten_Iten[0]=Ten_Iten[0]/norm(Ten_Iten[0]);
Ten_Iten[1]=Ten_Iten[1]/norm(Ten_Iten[1]);
Ten_Iten[2]=Ten_Iten[2]/norm(Ten_Iten[2]);
Ten_Iten[3]=Ten_Iten[3]/norm(Ten_Iten[3]);

return Ten_Iten;

}

auto Init_Tensors_one(vector<ITensor> & Ten_Iten, int& chi, int& D) -> vector<ITensor>
{
vector<ITensor> Ten_Iten1;
auto Ain=findtype(Ten_Iten[0],Atype);
auto Bout=findtype(Ten_Iten[0],Btype);

Ten_Iten1.push_back(randomTensor(Ain,prime(Ain,1),Bout,prime(Bout,1)));
Ten_Iten1[0]=Ten_Iten1[0]/norm(Ten_Iten1[0]);
Ten_Iten1.push_back(Ten_Iten[1]);
Ten_Iten1.push_back(Ten_Iten[2]);
Ten_Iten1.push_back(Ten_Iten[3]);


return Ten_Iten1;
}


auto Init_Tensors_fill(vector<ITensor> & Ten_Iten, int& chi, int& D, double & val) -> vector<ITensor>
{
for( auto i : range(Ten_Iten.size())  )
    Ten_Iten[i].fill(val);

return Ten_Iten;
}




auto Init_Corner(int& chi, int& D) -> vector<ITensor>
{
vector<ITensor> Env_Iten;


auto x=Index("x", chi, Xtype, 0);
auto y=Index("y", D, Ytype, 0);



Env_Iten.push_back(randomTensor(x,prime(x,1)));
Env_Iten.push_back(randomTensor(prime(x,1),y,prime(x,2)));
Env_Iten.push_back(randomTensor(prime(x,2),prime(y,1), prime(x,3)));
Env_Iten.push_back(randomTensor(prime(x,3),prime(x,4)));
/*4*/
Env_Iten.push_back(randomTensor(x,prime(y,3),prime(x,11)));
/*5-6*/
Env_Iten.push_back(randomTensor(prime(x,4),prime(y,2),prime(x,5)));
Env_Iten.push_back(randomTensor(prime(x,11),prime(y,4),prime(x,10)));
/*7-8*/
Env_Iten.push_back(randomTensor(prime(x,5),prime(y,5),prime(x,6)));
Env_Iten.push_back(randomTensor(prime(x,10),prime(x,9)));
/*9-11*/
Env_Iten.push_back(randomTensor(prime(x,7),prime(x,6)));
Env_Iten.push_back(randomTensor(prime(x,9),prime(y,7),prime(x,8)));
Env_Iten.push_back(randomTensor(prime(x,8),prime(y,6),prime(x,7)));


return Env_Iten;

}


auto norm_CTM( vector<ITensor> & Env_Iten,  vector<ITensor> & Ten_Iten) -> double
{

//x=noprime(Env_Iten[0].findtype(Xtype));
//auto y=noprime(findtype(Env_Iten[0],Ytype));
//auto Ain=noprime(findtype(Ten_Iten[0],Atype));
//auto Bout=noprime(findtype(Ten_Iten[0],Btype));
auto Ten_Iten1=Label_TensCTM(Env_Iten, Ten_Iten);

//println(Ten_Iten1[0]);

auto E1=Env_Iten[0] * Env_Iten[1] *Env_Iten[4] * Ten_Iten1[0] ;
auto E2=Env_Iten[3] * Env_Iten[2] *Env_Iten[5] *Ten_Iten1[1];
auto E3=Env_Iten[9] * Env_Iten[7] *Env_Iten[11] *Ten_Iten1[3];
auto E4=Env_Iten[8] * Env_Iten[10] *Env_Iten[6] *Ten_Iten1[2];
auto E_final=E1*E2*E3*E4;
//PrintData(E_final);
auto i=E_final.real();
return i;
}






auto Label_TensCTM( vector<ITensor> & Env_Iten,  vector<ITensor> & Ten_Iten) -> vector<ITensor>
{
vector<ITensor> Ten_Iten_tempo;
Ten_Iten_tempo.push_back(Ten_Iten[0]);
Ten_Iten_tempo.push_back(Ten_Iten[1]);
Ten_Iten_tempo.push_back(Ten_Iten[2]);
Ten_Iten_tempo.push_back(Ten_Iten[3]);

//x=noprime(Env_Iten[0].findtype(Xtype));
auto y=noprime(findtype(Env_Iten[1],Ytype));
auto Ain=noprime(findtype(Ten_Iten[0],Atype));
auto Bout=noprime(findtype(Ten_Iten[0],Btype));
//println(y);
//println(Ain);
//println(Bout);


Ten_Iten_tempo[0]*=delta(Ain,prime(y,3));
Ten_Iten_tempo[0]*=delta(Bout,prime(y,9));
Ten_Iten_tempo[0]*=delta(prime(Ain,1),prime(y,8));
Ten_Iten_tempo[0]*=delta(prime(Bout,1),y);

Ten_Iten_tempo[1]*=delta(Ain,prime(y,9));
Ten_Iten_tempo[1]*=delta(Bout,prime(y,10));
Ten_Iten_tempo[1]*=delta(prime(Ain,1),prime(y,2));
Ten_Iten_tempo[1]*=delta(prime(Bout,1),prime(y,1));

Ten_Iten_tempo[2]*=delta(Ain,prime(y,4));
Ten_Iten_tempo[2]*=delta(Bout,prime(y,7));
Ten_Iten_tempo[2]*=delta(prime(Ain,1),prime(y,11));
Ten_Iten_tempo[2]*=delta(prime(Bout,1),prime(y,8));

Ten_Iten_tempo[3]*=delta(Ain,prime(y,11));
Ten_Iten_tempo[3]*=delta(Bout,prime(y,6));
Ten_Iten_tempo[3]*=delta(prime(Ain,1),prime(y,5));
Ten_Iten_tempo[3]*=delta(prime(Bout,1),prime(y,10));

return Ten_Iten_tempo;
}


auto   Transpose_itensor( ITensor & E, Index & i, Index & j  )  -> ITensor
{
auto i_b=prime(j,1);
auto j_b=prime(i,1);
auto E1=ITensor(i_b,j_b);
randomize(E1);
for(auto s1 : range1(i.m()))
for(auto s2 : range1(j.m()))
{
E1.set(i_b(s2),j_b(s1),E.real(j(s2),i(s1)));
}
E1=prime(E1,i_b,-1);
E1=prime(E1,j_b,-1);
return E1;
}


auto Left_move(vector<ITensor> & Env_Iten, vector<ITensor> & Ten_Iten, int & chi) -> vector<ITensor>
{
auto Ten_Iten1=Label_TensCTM(Env_Iten, Ten_Iten);
auto x=noprime(findtype(Env_Iten[0],Xtype));
auto y=noprime(findtype(Env_Iten[1],Ytype));
float Inf = numeric_limits<float>::infinity();
wall_clock timer;
double n_secs;

/*******************************************************/
//println(Env_Iten[0]); 
auto E1=Env_Iten[0]*Env_Iten[1];
auto C = combiner(x,y);
E1*=C;
auto i_ind = commonIndex(C,E1);
auto j_ind=prime(x,2);
auto E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

auto E2=Env_Iten[8]*Env_Iten[10];
auto C1 = combiner(prime(x,10),prime(y,7));
E2*=C1;
auto i1_ind = commonIndex(C1,E2);
auto j1_ind=prime(x,8);
auto E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));

/*****************************************************/

E2=swapPrime(dag(E2),0,1);


E1=0.50*(E1+E2);
E1*=C;
E1*=prime(C,1);
ITensor U(x,y),D;
timer.tic();
Spectrum spec=diagHermitian(E1,U,D,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
n_secs = timer.toc();
//cout << "took " << n_secs << " seconds" << endl;
//Print(  norm(E1-dag(U)*D*prime(U)) );
Print(spec);
auto U_tran=U;
C = combiner(x,y);
U_tran*=C;
i_ind = commonIndex(C,U_tran);
j_ind=commonIndex(U_tran,D);
U_tran=Transpose_itensor(U_tran,i_ind,j_ind );
U_tran*=C;
//PrintData(U*prime(U_tran,j_ind,1));
U*=delta(x,prime(x,10));
U*=delta(y,prime(y,7));
/***************************************************************************/
/*******************************************************/

E1=Env_Iten[0]*Env_Iten[1]*Env_Iten[4]*Ten_Iten1[0];
C=combiner(prime(x,11),prime(y,8));
//println(C);
C1=combiner(prime(x,2),prime(y,9));
E1*=C;
E1*=C1;
i_ind = commonIndex(C,E1);
j_ind=commonIndex(C1,E1);
E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

E2=Env_Iten[8]*Env_Iten[10]*Env_Iten[6]*Ten_Iten1[2];
auto Com = combiner(prime(x,11),prime(y,8));
auto Com1 = combiner(prime(x,8),prime(y,11));
E2*=Com;
E2*=Com1;
i1_ind=commonIndex(Com,E2);
j1_ind=commonIndex(Com1,E2);
E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));
/*****************************************************/
E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1=E1*C;
E1=E1*prime(C,1);
ITensor U1(prime(x,11),prime(y,8)),D1;
spec=diagHermitian(E1,U1,D1,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);

auto U1_tran=U1;
C = combiner(prime(x,11),prime(y,8));
U1_tran*=C;
i_ind = commonIndex(C,U1_tran);
j_ind=commonIndex(U1_tran,D1);
U1_tran=Transpose_itensor(U1_tran,i_ind,j_ind );
U1_tran=U1_tran*C;
//PrintData(U1*prime(U1_tran,j_ind,1));
//PrintData(U1);
//PrintData(U1_tran);
/***************************************************************************/

vector<ITensor>  Env_Iten1=Env_Iten;

/**************************************************************************/
E1=Env_Iten[0]*Env_Iten[1]*U_tran;
auto d_link=findtype(E1,Link);
E1*=delta(d_link,x);
E1*=delta(prime(x,1),prime(x,2));
Env_Iten1[0]=E1/norm(E1);
E1=Env_Iten[8]*Env_Iten[10]*U;
d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x,10));
E1*=delta(prime(x,9),prime(x,8));
Env_Iten1[8]=E1/norm(E1);
/**************************************************************************/
U*=delta(prime(x,10),x);
U*=delta(prime(y,7),y);
auto dp_link=findtype(U,Link);
U=prime(U,dp_link,1);
E1=Env_Iten[4]*Ten_Iten1[0]*U1_tran*U;
E1*=delta(prime(y,9),prime(y,3));
d_link=findtype(U1_tran,Link);
E1*=delta(d_link,prime(x,11));
E1*=delta(prime(dp_link),x);
Env_Iten1[4]=E1/norm(E1);
//PrintData(E1);

U_tran*=delta(prime(x,10),x);
U_tran*=delta(prime(y,7),y);
dp_link=findtype(U1,Link);
U1=prime(U1,dp_link,1);
E1=Env_Iten[6]*Ten_Iten1[2]*U1*U_tran;
E1*=delta(prime(y,11),prime(y,4));
d_link=findtype(U_tran,Link);
E1*=delta(d_link,prime(x,10));
E1*=delta(prime(dp_link),prime(x,11));
Env_Iten1[6]=E1/norm(E1);
//PrintData(E1);
/**************************************************************************************/















/*******************************************************/
//println(Env_Iten[0]); 

E1=Env_Iten[2]*Env_Iten[3];
C = combiner(prime(x,4),prime(y,1));
E1*=C;
i_ind = commonIndex(C,E1);
j_ind=prime(x,2);
E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

E2=Env_Iten[9]*Env_Iten[11];
C1 = combiner(prime(x,6),prime(y,6));
E2*=C1;
i1_ind = commonIndex(C1,E2);
j1_ind=prime(x,8);
E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));

/*****************************************************/

E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1*=C;
E1*=prime(C,1);
//println(E1);
U=randomTensor(prime(x,4),prime(y,1));
timer.tic();
spec=diagHermitian(E1,U,D,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
n_secs = timer.toc();
//cout << "took " << n_secs << " seconds" << endl;

//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);
U_tran=U;
C = combiner(prime(x,4),prime(y,1));
U_tran*=C;
i_ind = commonIndex(C,U_tran);
j_ind=commonIndex(U_tran,D);
U_tran=Transpose_itensor(U_tran,i_ind,j_ind );
U_tran*=C;
//PrintData(U*prime(U_tran,j_ind,1));
U*=delta(prime(x,4),prime(x,6));
U*=delta(prime(y,1),prime(y,6));
/***************************************************************************/
/*******************************************************/

E1=Env_Iten[2]*Env_Iten[3]*Env_Iten[5]*Ten_Iten1[1];
C=combiner(prime(x,5),prime(y,10));
//println(C);
C1=combiner(prime(x,2),prime(y,9));
E1*=C;
E1*=C1;
i_ind = commonIndex(C,E1);
j_ind=commonIndex(C1,E1);
E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

E2=Env_Iten[7]*Env_Iten[9]*Env_Iten[11]*Ten_Iten1[3];
Com = combiner(prime(x,5),prime(y,10));
Com1 = combiner(prime(x,8),prime(y,11));
E2*=Com;
E2*=Com1;
i1_ind=commonIndex(Com,E2);
j1_ind=commonIndex(Com1,E2);
E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));
/*****************************************************/
E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1=E1*C;
E1=E1*prime(C,1);

U1=randomTensor(prime(x,5),prime(y,10));

spec=diagHermitian(E1,U1,D1,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);

U1_tran=U1;
C = combiner(prime(x,5),prime(y,10));
U1_tran*=C;
i_ind = commonIndex(C,U1_tran);
j_ind=commonIndex(U1_tran,D1);
U1_tran=Transpose_itensor(U1_tran,i_ind,j_ind );
U1_tran=U1_tran*C;
//PrintData(U1*prime(U1_tran,j_ind,1));
//PrintData(U1);
//PrintData(U1_tran);
/***************************************************************************/


/**************************************************************************/
E1=Env_Iten[2]*Env_Iten[3]*U_tran;
d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x,4));
E1*=delta(prime(x,3),prime(x,2));
Env_Iten1[3]=E1/norm(E1);



E1=Env_Iten[9]*Env_Iten[11]*U;
d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x,6));
E1*=delta(prime(x,7),prime(x,8));
Env_Iten1[9]=E1/norm(E1);



/**************************************************************************/
U*=delta(prime(x,6),prime(x,4));
U*=delta(prime(y,6),prime(y,1));
dp_link=findtype(U,Link);
U=prime(U,dp_link,1);
E1=Env_Iten[5]*Ten_Iten1[1]*U1_tran*U;
E1*=delta(prime(y,9),prime(y,2));
d_link=findtype(U1_tran,Link);
E1*=delta(d_link,prime(x,5));
E1*=delta(prime(dp_link),prime(x,4));
Env_Iten1[5]=E1/norm(E1);
//PrintData(E1);

U_tran*=delta(prime(x,4),prime(x,6));
U_tran*=delta(prime(y,1),prime(y,6));
dp_link=findtype(U1,Link);
U1=prime(U1,dp_link,1);
E1=Env_Iten[7]*Ten_Iten1[3]*U1*U_tran;
E1*=delta(prime(y,11),prime(y,5));
d_link=findtype(U_tran,Link);
E1*=delta(d_link,prime(x,6));
E1*=delta(prime(dp_link),prime(x,5));
Env_Iten1[7]=E1/norm(E1);
//PrintData(E1);
/**************************************************************************************/

return Env_Iten1;

}


auto Right_move(vector<ITensor> & Env_Iten, vector<ITensor> & Ten_Iten, int & chi) -> vector<ITensor>
{
auto Ten_Iten1=Label_TensCTM(Env_Iten, Ten_Iten);
auto x=noprime(findtype(Env_Iten[0],Xtype));
auto y=noprime(findtype(Env_Iten[1],Ytype));
float Inf = numeric_limits<float>::infinity();
wall_clock timer;
double n_secs;

/*******************************************************/
//println(Env_Iten[0]); 
auto E1=Env_Iten[0]*Env_Iten[4];
auto C = combiner(prime(x),prime(y,3));
E1*=C;
auto i_ind = commonIndex(C,E1);
auto j_ind=prime(x,11);
auto E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

auto E2=Env_Iten[3]*Env_Iten[5];
auto C1 = combiner(prime(x,3),prime(y,2));
E2*=C1;
auto i1_ind = commonIndex(C1,E2);
auto j1_ind=prime(x,5);
auto E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));

/*****************************************************/

E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1*=C;
E1*=prime(C,1);
ITensor U(prime(x),prime(y,3)),D;
timer.tic();
Spectrum spec=diagHermitian(E1,U,D,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
n_secs = timer.toc();
//cout << "took " << n_secs << " seconds" << endl;
//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);
auto U_tran=U;
C = combiner(prime(x),prime(y,3));
U_tran*=C;
i_ind = commonIndex(C,U_tran);
j_ind=commonIndex(U_tran,D);
U_tran=Transpose_itensor(U_tran,i_ind,j_ind );
U_tran*=C;
//PrintData(U*prime(U_tran,j_ind,1));
U*=delta(prime(x),prime(x,3));
U*=delta(prime(y,3),prime(y,2));
/***************************************************************************/
/*******************************************************/

E1=Env_Iten[0]*Env_Iten[1]*Env_Iten[4]*Ten_Iten1[0];
C=combiner(prime(x,2),prime(y,9));
//println(C);
C1=combiner(prime(x,11),prime(y,8));
E1*=C;
E1*=C1;
i_ind = commonIndex(C,E1);
j_ind=commonIndex(C1,E1);
E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

E2=Env_Iten[3]*Env_Iten[5]*Env_Iten[2]*Ten_Iten1[1];
auto Com = combiner(prime(x,2),prime(y,9));
auto Com1 = combiner(prime(x,5),prime(y,10));
E2*=Com;
E2*=Com1;
i1_ind=commonIndex(Com,E2);
j1_ind=commonIndex(Com1,E2);
E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));
/*****************************************************/
E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1=E1*C;
E1=E1*prime(C,1);
ITensor U1(prime(x,2),prime(y,9)),D1;
spec=diagHermitian(E1,U1,D1,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);

auto U1_tran=U1;
C = combiner(prime(x,2),prime(y,9));
U1_tran*=C;
i_ind = commonIndex(C,U1_tran);
j_ind=commonIndex(U1_tran,D1);
U1_tran=Transpose_itensor(U1_tran,i_ind,j_ind );
U1_tran=U1_tran*C;
//PrintData(U1*prime(U1_tran,j_ind,1));
//PrintData(U1);
//PrintData(U1_tran);
/***************************************************************************/

vector<ITensor>  Env_Iten1=Env_Iten;

/**************************************************************************/

E1=Env_Iten[0]*Env_Iten[4]*U_tran;
auto d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x));
E1*=delta(prime(x,11),x);
Env_Iten1[0]=E1/norm(E1);


E1=Env_Iten[3]*Env_Iten[5]*U;
d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x,3));
E1*=delta(prime(x,5),prime(x,4));
Env_Iten1[3]=E1/norm(E1);
/**************************************************************************/
U*=delta(prime(x,1),prime(x,3));
U*=delta(prime(y,3),prime(y,2));
auto dp_link=findtype(U,Link);
U=prime(U,dp_link,1);
E1=Env_Iten[1]*Ten_Iten1[0]*U1_tran*U;
E1*=delta(prime(y,8),y);
d_link=findtype(U1_tran,Link);
E1*=delta(d_link,prime(x,2));
E1*=delta(prime(dp_link),prime(x));
Env_Iten1[1]=E1/norm(E1);
//PrintData(E1);



U_tran*=delta(prime(x,3),prime(x));
U_tran*=delta(prime(y,2),prime(y,3));
dp_link=findtype(U1,Link);
U1=prime(U1,dp_link,1);
E1=Env_Iten[2]*Ten_Iten1[1]*U1*U_tran;
E1*=delta(prime(y,10),prime(y,1));
d_link=findtype(U_tran,Link);
E1*=delta(d_link,prime(x,3));
E1*=delta(prime(dp_link),prime(x,2));
Env_Iten1[2]=E1/norm(E1);


//PrintData(E1);
/**************************************************************************************/













/*******************************************************/
//println(Env_Iten[0]); 

E1=Env_Iten[8]*Env_Iten[6];
C = combiner(prime(x,9),prime(y,4));
E1*=C;
i_ind = commonIndex(C,E1);
j_ind=prime(x,11);
E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

E2=Env_Iten[9]*Env_Iten[7];
C1 = combiner(prime(x,7),prime(y,5));
E2*=C1;
i1_ind = commonIndex(C1,E2);
j1_ind=prime(x,5);
E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));

/*****************************************************/

E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1*=C;
E1*=prime(C,1);
//println(E1);
U=randomTensor(prime(x,9),prime(y,4));
timer.tic();
spec=diagHermitian(E1,U,D,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
n_secs = timer.toc();
//cout << "took " << n_secs << " seconds" << endl;

//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);
U_tran=U;
C = combiner(prime(x,9),prime(y,4));
U_tran*=C;
i_ind = commonIndex(C,U_tran);
j_ind=commonIndex(U_tran,D);
U_tran=Transpose_itensor(U_tran,i_ind,j_ind );
U_tran*=C;
//PrintData(U*prime(U_tran,j_ind,1));
U*=delta(prime(x,9),prime(x,7));
U*=delta(prime(y,4),prime(y,5));
/***************************************************************************/
/*******************************************************/

E1=Env_Iten[8]*Env_Iten[6]*Env_Iten[10]*Ten_Iten1[2];
C=combiner(prime(x,8),prime(y,11));
//println(C);
C1=combiner(prime(x,11),prime(y,8));
E1*=C;
E1*=C1;
i_ind = commonIndex(C,E1);
j_ind=commonIndex(C1,E1);
E1_trans=Transpose_itensor(E1,i_ind,j_ind );
E1*=prime(E1_trans,i_ind,1);

/*****************************************************/

E2=Env_Iten[7]*Env_Iten[9]*Env_Iten[11]*Ten_Iten1[3];
Com = combiner(prime(x,8),prime(y,11));
Com1 = combiner(prime(x,5),prime(y,10));
E2*=Com;
E2*=Com1;
i1_ind=commonIndex(Com,E2);
j1_ind=commonIndex(Com1,E2);
E2_trans=Transpose_itensor(E2,i1_ind,j1_ind );
E2*=prime(E2_trans,i1_ind,1);
E2*=delta(i1_ind, i_ind);
E2*=delta(prime(i1_ind), prime(i_ind));
/*****************************************************/
E2=swapPrime(dag(E2),0,1);
E1=0.50*(E1+E2);
E1=E1*C;
E1=E1*prime(C,1);

U1=randomTensor(prime(x,8),prime(y,11));

spec=diagHermitian(E1,U1,D1,{"Cutoff=",Inf,"Minm=",chi,"IndexType=",Link});
//Print(  norm(E1-dag(U)*D*prime(U)) );
//Print(spec);

U1_tran=U1;
C = combiner(prime(x,8),prime(y,11));
U1_tran*=C;
i_ind = commonIndex(C,U1_tran);
j_ind=commonIndex(U1_tran,D1);
U1_tran=Transpose_itensor(U1_tran,i_ind,j_ind );
U1_tran=U1_tran*C;
//PrintData(U1*prime(U1_tran,j_ind,1));
//PrintData(U1);
//PrintData(U1_tran);
/***************************************************************************/


/**************************************************************************/
E1=Env_Iten[8]*Env_Iten[6]*U_tran;
d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x,9));
E1*=delta(prime(x,11),prime(x,10));
Env_Iten1[8]=E1/norm(E1);



E1=Env_Iten[9]*Env_Iten[7]*U;
d_link=findtype(E1,Link);
E1*=delta(d_link,prime(x,7));
E1*=delta(prime(x,6),prime(x,5));
Env_Iten1[9]=E1/norm(E1);



/**************************************************************************/
U*=delta(prime(x,7),prime(x,9));
U*=delta(prime(y,5),prime(y,4));
dp_link=findtype(U,Link);
U=prime(U,dp_link,1);
E1=Env_Iten[10]*Ten_Iten1[2]*U1_tran*U;
E1*=delta(prime(y,8),prime(y,7));
d_link=findtype(U1_tran,Link);
E1*=delta(d_link,prime(x,8));
E1*=delta(prime(dp_link),prime(x,9));
Env_Iten1[10]=E1/norm(E1);

U_tran*=delta(prime(x,7),prime(x,9));
U_tran*=delta(prime(y,5),prime(y,4));
dp_link=findtype(U1,Link);
U1=prime(U1,dp_link,1);
E1=Env_Iten[11]*Ten_Iten1[3]*U1*U_tran;
E1*=delta(prime(y,6),prime(y,10));
d_link=findtype(U_tran,Link);
E1*=delta(d_link,prime(x,7));
E1*=delta(prime(dp_link),prime(x,8));
Env_Iten1[11]=E1/norm(E1);
/**************************************************************************************/
return Env_Iten1;
}


auto Permute_Env(vector<ITensor> & Env_left) -> vector<ITensor>
{
auto x=noprime(findtype(Env_left[0],Xtype));
auto y=noprime(findtype(Env_left[1],Ytype));

vector<ITensor> Env_left1;
Env_left1=Env_left;


Env_left1[2]=Env_left[1];
Env_left1[2]*=delta(prime(x,3),prime(x,2));
Env_left1[2]*=delta(prime(x,2),prime(x,1));
Env_left1[2]*=delta(prime(y,1),y);

Env_left1[1]=Env_left[2];
Env_left1[1]*=delta(prime(x,2),prime(x,1));
Env_left1[1]*=delta(prime(x,3),prime(x,2));
Env_left1[1]*=delta(prime(y,1),y);


Env_left1[11]=Env_left[10];
Env_left1[11]*=delta(prime(x,8),prime(x,7));
Env_left1[11]*=delta(prime(x,9),prime(x,8));
Env_left1[11]*=delta(prime(y,6),prime(y,7));

Env_left1[10]=Env_left[11];
Env_left1[10]*=delta(prime(x,9),prime(x,8));
Env_left1[10]*=delta(prime(x,8),prime(x,7));
Env_left1[10]*=delta(prime(y,6),prime(y,7));

return Env_left1;
}
auto Permute_Ten(vector<ITensor> & Ten_Iten) -> vector<ITensor>
{
auto Ain=findtype(Ten_Iten[0],Atype);
auto Bout=findtype(Ten_Iten[0],Btype);

vector<ITensor> Ten_Iten1;
Ten_Iten1=Ten_Iten;

Ten_Iten1[0]=Ten_Iten[1];
Ten_Iten1[1]=Ten_Iten[0];

Ten_Iten1[2]=Ten_Iten[3];
Ten_Iten1[3]=Ten_Iten[2];


return Ten_Iten1;
}


auto Permute1_Env(vector<ITensor> & Env_left) -> vector<ITensor>
{
auto x=noprime(findtype(Env_left[0],Xtype));
auto y=noprime(findtype(Env_left[1],Ytype));

vector<ITensor> Env_left1;
Env_left1=Env_left;

Env_left1[6]=Env_left[4];
Env_left1[6]*=delta(prime(x,10),prime(x,11));
Env_left1[6]*=delta(prime(x,11),x);
Env_left1[6]*=delta(prime(y,4),prime(y,3));

Env_left1[4]=Env_left[6];
Env_left1[4]*=delta(prime(x,11),x);
Env_left1[4]*=delta(prime(x,10),prime(x,11));
Env_left1[4]*=delta(prime(y,4),prime(y,3));


Env_left1[5]=Env_left[7];
Env_left1[5]*=delta(prime(x,5),prime(x,4));
Env_left1[5]*=delta(prime(x,6),prime(x,5));
Env_left1[5]*=delta(prime(y,2),prime(y,5));

Env_left1[7]=Env_left[5];
Env_left1[7]*=delta(prime(x,6),prime(x,5));
Env_left1[7]*=delta(prime(x,5),prime(x,4));
Env_left1[7]*=delta(prime(y,2),prime(y,5));

return Env_left1;
}
auto Permute1_Ten(vector<ITensor> & Ten_Iten) -> vector<ITensor>
{
auto Ain=findtype(Ten_Iten[0],Atype);
auto Bout=findtype(Ten_Iten[0],Btype);

vector<ITensor> Ten_Iten1;
Ten_Iten1=Ten_Iten;

Ten_Iten1[0]=Ten_Iten[2];
Ten_Iten1[2]=Ten_Iten[0];

Ten_Iten1[1]=Ten_Iten[3];
Ten_Iten1[3]=Ten_Iten[1];


return Ten_Iten1;
}



