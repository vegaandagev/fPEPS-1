set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"
set key t r 
set key box 6
set key font ",40"

set key spacing 6
set output "Corr.eps"


set xtics ( 0.5, 1, 1.4);
set ytics ( -1.0," " -2, -3," " -4, -5," " -6,-7);


p  [0.0:10]   "CorrelationHH5.txt"   t 'D=5'  with points  pointtype 2 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationHH6.txt"   t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371", "CorrelationHH7.txt" t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED", "CorrelationVV5.txt" t 'D=5'  with points  pointtype 2 lw 3  ps 5.0 lc rgb "red","CorrelationVV6.txt" t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "red", "CorrelationVV7.txt"  t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "red"



