set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel "C^{d}"   font ",40" textcolor rgb "black"

set key b l 
set key box 6
set key font ",45"
set key spacing 6
set output "Corrlog.eps"

set xrange [log(10)/log(10):log(25)/log(10)]

set xtics (  1,1.4,3, 4, 5);
set ytics nomirror ( -3, -5, -6, -14);



#set ytics nomirror autofreq
#set border linewidth 1

#set format ‘$%g$’
#set format y ‘$10^{%T}$’
#set logscale x
#set logscale y
M(x)=(c1*x)+d1

fit  [log(1)/log(10):log(4)/log(10)] M(x) "CorrelationVV7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    via c1,d1 

p  [0.01:1.6] "CorrelationVV5.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t "D=5"   with points  pointtype 2 lw 3  ps 5.0 lc rgb "#4B0082", "CorrelationVV6.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=6'   with points  pointtype 4 lw 3  ps 5.0 lc rgb "#800000", "CorrelationVV7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#191970", "CorrelationINF.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D {/Symbol \256 } {/Symbol \245 }' with points  pointtype 6 lw 3  ps 5.0 lc rgb "#4682B4",  M(x) t 'C(r)=r^{-2.9}' with line dt 4 lw 4 linetype 4 lc rgb "blue"



#"CorrelationHF.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t 'D {/#Symbol \256 } {/Symbol \245 }' pointtype 11  lc rgb "blue" ps #3, 
#"CorrelationVV8.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) t "D=8" with points  pointtype 12 lw 3  ps 5.0 lc rgb "#FF6347"
