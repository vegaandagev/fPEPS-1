set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/D"   font ",40" textcolor rgb "black"
set ylabel "m, {/Symbol D}T"  font ",40"  textcolor rgb "black"

set xtics nomirror ("1/5" .2, "1/9" .111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.473,-0.477,-0.4782, -0.482);
set tics scale 3


set key box c l  
set key font ",40"
set key spacing 8
set key width 1

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

#G(x)=a+b*(x**(2))+c*(x**(4))
#fit [0.10:0.33]  G(x)  "ipeps.txt" u ($1):($3)    via a,b,c
M(x)=a1+b1*(x)
fit [0.10:0.33]  M(x)  "ipeps.txt" u ($1):($3)    via a1,b1#,c
#M1(x)=a2+b2*(x**(1))+c2*(x**(2))
#fit   M1(x) "J20.5.txt" u ($1):($2)    via a2,b2,c2
#P(x)=a1+b1*(x**(2))+c1*(x**(4))
#fit   P(x) "DMRG.txt" u ($1):($2)    via a1,b1,c1
#K(x)=a3+b3*(x**(2))+c3*(x**(4))
#fit   K(x) "ipeps.txt" u ($1):($2)    via a3,b3,c3



G1(x)=0


set output "Mag.eps"
p [0.0:0.40]   "ipeps.txt" u ($1):($3) t "m" with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B","ipeps.txt" u ($1):($5) t "{/Symbol D}T_{x}" with points  pointtype 4 lw 3  ps 8.0 lc rgb "#C71585","ipeps.txt" u ($1):($6) t "{/Symbol D}T_{x-y}" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#FF4500", x<0.330 ? M(x) : 1/0 notitle lw 4 lc rgb "#DC143C" dt 4


