set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40”
set ytics font "Times-Roman, 40”

set xlabel "1/D"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}T_{y}"  font ",40"  textcolor rgb "black"
#set log y 
set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
set ytics (0.4,0.2,0.3,0.1,0.05,0.01)
set tics scale 4
set xtics nomirror


Z1(x)=a4+b4*(x**(1))
fit   Z1(x) "J20.55.txt" u ($1):($2)    via a4,b4#,c3

Z2(x)=a5+b5*(x**(1))
fit   Z2(x) "J20.6.txt" u ($1):($2)    via a5,b5#,c5

Z3(x)=a1+b1*(x**(1))
fit   Z3(x) "J20.54.txt" u ($1):($2)    via a1,b1#,c5

set key box 2
set key t l 
set key font ",30"
set key spacing 6

#set key at 0.1,0.7
set output "TI.eps"
set multiplot
p  [0.008:0.30] [0.0:0.36] 'J20.6.txt'  t 'J2=0.60'  with points  pointtype 3 lw 3  ps 6.0 lc rgb "#FF1493",'J20.55.txt'  t 'J2=0.55'   with points  pointtype 2 lw 3  ps 6.0 lc rgb "#8B008B", 'J20.54.txt'  t 'J2=0.54'  with points  pointtype 4 lw 3  ps 6.0 lc rgb "#8A2BE2", a4+b4*(x**(1)) notitle lw 4 lc rgb "#8B008B" dt 4, a5+b5*(x**(1)) notitle lw 4 lc rgb "#8A2BE2" dt 4, a1+b1*(x**(1)) notitle lw 4 lc rgb "#8A2BE2" dt 4

set origin .28, .15
set size 0.55,0.320
clear
unset key
unset xlabel
set ylabel "{/Symbol D}T_{x}"  font ",30"  textcolor rgb "black"
set ytics font "Times-Roman, 30”
set xtics font "Times-Roman, 30”
set tics scale 2
set xtics nomirror

set ytics (0.4,0.2,0.3,0.1,0.01,0.05)


p  [0.04:0.34] [0.0:0.1]  'J20.6.txt' u ($1):($3)  t 'J2=0.60'  with points  pointtype 3 lw 2  ps 4.0 lc rgb "#FF1493",'J20.55.txt' u ($1):($3)  t 'J2=0.55'   with points  pointtype 2 lw 2  ps 4.0 lc rgb "#8B008B", 'J20.54.txt' u ($1):($3)  t 'J2=0.54'  with points  pointtype 4 lw 2  ps 4.0 lc rgb "#8A2BE2"
unset multiplot

