set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"
set key b l 
set key box 6
set key font ",40"



set key spacing 6


set output "Corrlog.eps"

#set xrange [log(2)/log(10):log(45)/log(10)]

set xtics ( 0.5, 1, 1.4);
set ytics ( -1.0," " -2, -3," " -4, -5," " -6,-7);

M(x)=(c1*x)+d1

fit  [log(1)/log(10):log(5)/log(10)] M(x) "CorrelationH7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    via c1,d1 

p  [0.0:1.4]   "CorrelationH5.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=5'  with points  pointtype 2 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationH6.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371", "CorrelationH7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED", "CorrelationV5.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=5'  with points  pointtype 2 lw 3  ps 5.0 lc rgb "red","CorrelationV6.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "red", "CorrelationV7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "red",  M(x) t 'C(r)=r^{-1.65}' with line dt 4 lw 4 linetype 4 lc rgb "blue"

#, "CorrelationH8.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t #'D=8' with points  pointtype 12 lw 3  ps 5.0 lc rgb "#FF6347", "CorrelationHF.txt"  u #((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D {/Symbol \256 } {/Symbol #\245 }' with points  pointtype 14 lw 3  ps 5.0 lc rgb "#9400D3",


