#include "Header.h"
int main()
{

//println("hi");
int chi=20;
int D=2;
wall_clock timer;
double n_secs;
double val=1.00;
double norm;
double norm1;
srand ( 0 );
auto Env_Iten=Init_Corner(chi, D);
auto Ten_Iten=Init_Tensors(chi, D);
auto Ten_Iten1=Init_Tensors_one(Ten_Iten,chi, D);
//Ten_Iten=Init_Tensors_fill(Ten_Iten,chi, D,val);


//PrintData(Ten_Iten[0]); 
//PrintData(Ten_Iten1[0]); 

//timer.tic();
//auto norm=norm_CTM(Env_Iten, Ten_Iten);
//n_secs = timer.toc();
//cout << "took " << n_secs << " seconds" << endl;
//auto norm1=norm_CTM(Env_Iten, Ten_Iten1);
cout<<setprecision(12)<<endl;
//println(norm1);

//Env_Iten=Init_Tensors_fill(Env_Iten,chi, D,val);



auto Env_left =Env_Iten;

for (auto i : range1(10))
{
    println(i);
    Env_left =Left_move(Env_left,Ten_Iten,chi);
    Env_left=Permute_Env(Env_left);
    Ten_Iten=Permute_Ten(Ten_Iten);
    Env_left =Left_move(Env_left,Ten_Iten,chi);
    Env_left=Permute_Env(Env_left);
    Ten_Iten=Permute_Ten(Ten_Iten);

    Env_left =Right_move(Env_left,Ten_Iten,chi);
    Env_left=Permute1_Env(Env_left);
    Ten_Iten=Permute1_Ten(Ten_Iten);

    Env_left =Right_move(Env_left,Ten_Iten,chi);
    Env_left=Permute1_Env(Env_left);
    Ten_Iten=Permute1_Ten(Ten_Iten);

    norm=norm_CTM(Env_left, Ten_Iten);
    norm1=norm_CTM(Env_left, Ten_Iten1);
    println(norm/norm1);

    //PrintData(Env_left[0]);

}


    return 0;
    }


