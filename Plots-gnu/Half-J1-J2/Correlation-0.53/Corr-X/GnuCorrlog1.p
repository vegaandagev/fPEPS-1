#set terminal postscript eps enhanced color font 'Helvetica,14' lw 3
set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "log(r)"   font ",40" textcolor rgb "black"
set ylabel 'log(C(r))'  font ",40"  textcolor rgb "black"
set key b l 
set key box 6
set key font ",45"

#set key width 6

set key spacing 6


set output "Corrlog1.eps"

#set xrange [log(2)/log(10):log(45)/log(10)]

set xtics ( 0.5, 1, 1.4);
set ytics ( -1.0, -5, -9, -12);


p       "CorrelationH5.txt"  u ($1):((log(abs($2)))*(1/log(10))) t 'D=5'  with points  pointtype 2 lw 3  ps 5.0 lc rgb "#8B008B", "CorrelationH6.txt"  u ($1):((log(abs($2)))*(1/log(10))) t 'D=6'  with points  pointtype 3 lw 3  ps 5.0 lc rgb "#800000", "CorrelationH7.txt"  u ($1):((log(abs($2)))*(1/log(10))) t 'D=7' with points  pointtype 4 lw 3  ps 5.0 lc rgb "#191970", "CorrelationH8.txt"  u ($1):((log(abs($2)))*(1/log(10))) t 'D=8' with points  pointtype 9 lw 3  ps 5.0 lc rgb "#191970"



