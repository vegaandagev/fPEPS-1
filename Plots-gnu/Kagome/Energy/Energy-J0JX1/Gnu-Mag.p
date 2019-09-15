set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/D"   font ",35" textcolor rgb "black"
set ylabel "M"  font ",35"  textcolor rgb "black"

set xtics nomirror ("1/5" .2,"1/12" 0.083,"1/8" 0.125)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key box t l  
set key font ",35"
set key spacing 6

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

#G(x)=a+b*(x**(2))+c*(x**(4))
#fit G(x) "Echiral5.txt" u ($1):($3)    via a,b,c

M(x)=b1*((x**c1))
fit  M(x) "Echiral.txt" u ($1):($3)    via b1,c1

G(x)=b*((x**c))
fit  G(x) "EchiralS.txt" u ($1):($3)    via b,c


#[0.0:0.16]
#P(x)=exp(-x*d4)
#fit  P(x) "Echiral5.txt" u ($1):($3)    via d4


#

set output "M.eps"
p [0.0:0.250] [0:0.15]     "Echiral.txt"  u ($1):($3)  t "J=0:uniform"  with points  pointtype 2 lw 3  ps 8.0 lc rgb "#cc0000","EchiralS.txt" u  ($1):($3) t "J=0:Staggered" with points  pointtype 4 lw 3  ps 8.0 lc rgb "#8b008b", x<0.27 ? M(x) : 1/0 notitle lw 4 lc rgb "#cc0000" dt 4, x<0.335 ? G(x) : 1/0 notitle lw 4 lc rgb "#8b008b" dt 4

#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,

