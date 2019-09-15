set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/D"   font ",40" textcolor rgb "black"
set ylabel "m"  font ",40"  textcolor rgb "black"

set xtics ("1/5" .2, "1/8" .12511,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.473,-0.477,-0.4782, -0.482);
set tics scale 3


set key box t l  
set key font ",40"
set key spacing 6

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(1))
fit [0.10:0.31]  G(x) "55.txt" u ($1):($3)    via a,b
M(x)=a1+b1*(x**(1))
fit [0.10:0.31]  M(x) "545.txt" u ($1):($3)    via a1,b1
#C(x)=a2+b2*(x**(1))
#fit [0.10:0.31]  C(x) "548.txt" u ($1):($3)    via a2,b2

#M1(x)=a2+b2*(x)
#fit [0.10:0.31]  M1(x) "ipeps.txt" u ($1):($5)    via a2,b2
#P(x)=a1+b1*(x**(2))+c1*(x**(4))
#fit  [0.10:0.34] P(x) "ipeps.txt" u ($1):($4)    via a1,b1,c1
#K(x)=a3+b3*(x**(2))+c3*(x**(4))
#fit   K(x) "ipeps.txt" u ($1):($2)    via a3,b3,c3


#G1(x)=0
#M(x)=0
#C(x)=0

set output "Mag.eps"
p [0.0:0.40] [0.35:1.0]   x<0.330 ? M(x) : 1/0 notitle lw 4 lc rgb "#DC143C" dt 4,x<0.330 ? G(x) : 1/0 notitle lw 4 lc rgb "#FF00FF" dt 4,"55.txt" u ($1):($2)  t 'Stripe:J_2=0.55' with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B","545.txt" u ($1):($3)  t 'Stripe:J_2=0.545' with points  pointtype 8 lw 3  ps 8.0 lc rgb "#1E90FF","55.txt" u ($1):($3)  t 'Neel:J_2=0.55' with points  pointtype 4 lw 3  ps 8.0 lc rgb "#FF1493","545.txt" u ($1):($2)  t 'Neel:J_2=0.545' with points  pointtype 10 lw 3  ps 8.0 lc rgb "#8B0000"


#, x<0.255 ? P(x) : 1/0 notitle lw 4 lc rgb "#C71585" dt 4
#, "DMRG.txt"  t "DMRG" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#FF1493"
#"548.txt" u ($1):($2)  t 'Stripe:J_2=0.548' with points  pointtype 6 lw 3  ps 8.0 lc rgb "#6A5ACD","548.txt" u ($1):($2)  t 'Stripe:J_2=0.548' with points  pointtype 6 lw 3  ps 8.0 lc rgb "#6A5ACD"

