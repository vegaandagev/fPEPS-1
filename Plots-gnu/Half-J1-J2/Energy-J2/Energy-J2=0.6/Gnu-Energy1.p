set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/D or 1/L_y"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"

set xtics ("1/5" .2, "1/9" .111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
set ytics (-0.473,-0.477,-0.4782, -0.482);


set key box t l  
set key font ",40"
set key spacing 6
set tics scale 4
set xtics nomirror

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(2))+c*(x**(4))
fit [0.10:0.66]  G(x) "J20.5.txt" u ($1):($2)    via a,b,c
M(x)=a1+b1*(x)
fit [0.10:0.167]  M(x) "J20.5.txt" u ($1):($2)    via a1,b1#,c
#M1(x)=a2+b2*(x**(1))+c2*(x**(2))
#fit   M1(x) "J20.5.txt" u ($1):($2)    via a2,b2,c2
P(x)=a1+b1*(x**(2))#+c1*(x**(4))
fit   P(x) "CDMRG.txt" u ($1):($2)    via a1,b1#,c1
K(x)=a3+b3*(x**(2))+c3*(x**(4))
fit   K(x) "PEPS.txt" u ($1):($2)    via a3,b3,c3
 
#set x2label '1/L_y' font ",40" textcolor rgb "black"
#set x2tics ("1/5" .2, "1/8" .125,"1/3" .333)

G1(x)=0


set output "HeisenbergJ2-5.eps"
p [0.0:0.40] [-0.482:-0.4710]   -0.4773 t "Cluster PEPS" with line lw 4 linetype 6, x<0.2550 ? G(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4,   "J20.5.txt"  notitle with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B","CDMRG.txt"  t "DMRG" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#FF1493", x<0.255 ? P(x) : 1/0 notitle lw 4 lc rgb "#C71585" dt 4


