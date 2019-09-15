set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/D"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"
set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)

set key t r 
set key font ",30"
set key box 6
set key font ",40"
set tics scale 4
set key spacing 6



set output "Corrlog.eps"

M(x)=(c1*(x**(1)))+d1
fit  M(x) "fit.txt"  via c1,d1 
M1(x)=(c2*(x**(1)))+d2
fit  M1(x) "fit1.txt"  via c2,d2 
M2(x)=(c3*(x**(1)))+d3
fit  M2(x) "fit2.txt"  via c3,d3 
M3(x)=(c4*(x**(1)))+d4
fit  M3(x) "fit3.txt"  via c4,d4 

p  [0:0.4]  "fit1.txt" t "r=6"  with points  pointtype 4 lw 3  ps 8.0 lc rgb "#4682B4", "fit3.txt"  t "r=8"  with points  pointtype 6 lw 3  ps 8.0 lc rgb "#6495ED","fit2.txt" t "r=10"  with points  pointtype 8 lw 3  ps 8.0 lc rgb "#FF6347","fit.txt"  t "r=12"  with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B", x<0.27 ? M(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4, x<0.25 ? M1(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4, x<0.25 ? M2(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4,x<0.25 ? M3(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4




