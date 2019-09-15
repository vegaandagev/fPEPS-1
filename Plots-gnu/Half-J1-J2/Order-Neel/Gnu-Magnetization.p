set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40”
set ytics font "Times-Roman, 40”

set xlabel "1/D"   font ",40" textcolor rgb "black"
set ylabel "m"  font ",40"  textcolor rgb "black"

set xtics nomirror ("1/5" .2, "1/9" .111,"1/3" .333)
set ytics nomirror (0.05,0.2,0.3,0.1,0)
set tics scale 4

set key box 2

M(x)=a1+b1*(x**(1))#+c1*(x**(2))
fit   M(x) "J20.txt" u ($1):($2)    via a1,b1#,c1

Z(x)=a3+b3*(x**(1))#+c3*(x**(2))
fit  Z(x) "J20.45.txt" u ($1):($2)    via a3,b3

Z1(x)=a4+b4*(x**(1))
fit   Z1(x) "J20.55.txt" u ($1):($2)    via a4,b4

H(x)=a2+b2*(x**(1))
fit   H(x) "J20.50.txt" u ($1):($2)    via a2,b2

H1(x)=a6+b6*(x**(1))
fit  [0.01:0.24] H1(x) "J20.53.txt" u ($1):($2)    via a6,b6


set key t l 
set key font ",35”
set key spacing 6

set output "MagnetizationJ2-5.eps"

p  [0.0:0.35] [0.0:0.49]  'J20.txt'  t 'J_{2}=0.0' with points  pointtype 2 lw 3  ps 6.0 lc rgb "#8B008B",'J20.45.txt'  t 'J_{2}=0.45' with points  pointtype 3 lw 3  ps 6.0 lc rgb "#FF1493",'J20.50.txt'  t 'J_{2}=0.50'  with points  pointtype 4 lw 3  ps 6.0 lc rgb "#8A2BE2",'J20.53.txt'  t 'J_{2}=0.53'  with points  pointtype 6 lw 3  ps 6.0 lc rgb "#DC143C", a1+b1*(x**(1)) notitle lw 4 lc rgb "#8B008B" dt 4, a2+b2*(x**(1)) notitle lw 4 lc rgb "#8A2BE2" dt 4, a3+b3*(x**(1)) notitle lw 4 lc rgb "#FF1493" dt 4,  a6+b6*(x**(1)) notitle lw 4 lc rgb "#DC143C" dt 4


