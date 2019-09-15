#set terminal postscript eps enhanced color font 'Helvetica,14' lw 3
set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}(r)'  font ",40"  textcolor rgb "black"
set key b l 
set key box 6
set key font ",40"

#set key width 6

set key spacing 6


set output "Corrlog1.eps"

set xrange [log(2)/log(2):log(45)/log(2)]

set xtics ( 5, 15, 20,25,30,35, 40);
set ytics ( -3,-5,-7, -9,-11, -13,-15);


M4(x)=(c4*x)+d4
fit  [10:40] M4(x) "CorrelationHH4.txt" u ($1):((log(abs($2)))*(1/1))    via c4,d4 

M5(x)=(c5*x)+d5
fit  [10:40] M5(x) "CorrelationHH5.txt" u ($1):((log(abs($2)))*(1/1))    via c5,d5 

M6(x)=(c6*x)+d6
fit  [10:40] M6(x) "CorrelationVV6.txt" u ($1):((log(abs($2)))*(1/1))    via c6,d6 

M7(x)=(c7*x)+d7
fit  [10:40] M7(x) "CorrelationHH7.txt" u ($1):((log(abs($2)))*(1/1))    via c7,d7 


M8(x)=(c8*x)+d8
fit  [10:40] M8(x) "CorrelationVV8.txt" u ($1):((log(abs($2)))*(1/1))    via c8,d8 




p  [5:40]  "CorrelationHH4.txt"  u ($1):((log(abs($2)))*(1/1)) t 'D=4'  with points  pointtype 6 lw 3  ps 5.0 lc rgb "#48D1CC", "CorrelationHH5.txt"  u ($1):((log(abs($2)))*(1/1)) t 'D=5'  with points  pointtype 5 lw 3  ps 5.0 lc rgb "#FF69B4", "CorrelationVV6.txt"  u ($1):((log(abs($2)))*(1/1)) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationHH7.txt"  u ($1):((log(abs($2)))*(1/1)) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED", "CorrelationVV8.txt"  u ($1):((log(abs($2)))*(1/1)) t 'D=8' with points pointtype 12 lw 3  ps 5.0 lc rgb "#FF6347", M4(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", M5(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", M6(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",  M7(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",  M8(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970"



