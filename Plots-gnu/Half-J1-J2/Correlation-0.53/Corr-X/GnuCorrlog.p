#set terminal postscript eps enhanced color font 'Helvetica,14' lw 3
set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C(r)'  font ",40"  textcolor rgb "black"
set key b l 
set key box 6
set key font ",45"

#set key width 6

set key spacing 6


set output "Corrlog.eps"

set xrange [log(2)/log(10):log(45)/log(10)]

set xtics ( 0.5, 1, 1.4);
set ytics ( -1.0, -5, -9, -12);

M(x)=(c1*x)+d1

fit  [log(2)/log(10):log(7)/log(10)] M(x) "CorrelationH7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    via c1,d1 

p  [log(1)/log(10):log(30)/log(10)]   "CorrelationH5.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=5'  with points  pointtype 3 lw 3  ps 5.0 lc rgb "#8B008B", "CorrelationH6.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371", "CorrelationH7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#191970", M(x) t 'C(r)=r^{-2.0}' with line dt 4 lw 4 linetype 4 lc rgb "blue"

#, "CorrelationH7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED", #"CorrelationHF.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D {/Symbol \256 } {/Symbol \245 }' with points  pointtype 6 lw 3  #ps 5.0 lc rgb "#4682B4",


