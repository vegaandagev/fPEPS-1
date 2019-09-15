set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel "C^{d}"   font ",40" textcolor rgb "black"

set key b l 
set key box 6
set key font ",40”
set key spacing 6



set logscale x
set logscale y


set tics nomirror out scale 2.5


set format "%g"
set format y "10^{%T}"

#set ytics ( 10e-2, 10e-4, 10e-6, 10e-8);

set xtics nomirror autofreq font ",30"
set ytics nomirror autofreq font ",30"


M(x)=(c1*x)+d1

fit  [log(2)/log(2):log(5)/log(2)] M(x) "CorrelationVV7.txt" u ((log(abs($1)))*(1/log(2))):((log(abs($2)))*(1/log(2)))    via c1,d1
set output "Corr1.eps"
p  [1.9:20] "CorrelationVV5.txt"  u (((abs($1)))*(1/(1))):(((abs($2)))*(1/(1))) t "D=5"   with points  pointtype 2 lw 3  ps 5.0 lc rgb "#4B0082", "CorrelationVV6.txt"  u (((abs($1)))*(1/(1))):(((abs($2)))*(1/(1))) t 'D=6'   with points  pointtype 4 lw 3  ps 5.0 lc rgb "#800000", "CorrelationVV7.txt"  u (((abs($1)))*(1/(1))):(((abs($2)))*(1/(1))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#191970", "CorrelationVV8.txt"  u (((abs($1)))*(1/(1))):(((abs($2)))*(1/(1))) t "D=8" with points  pointtype 12 lw 3  ps 5.0 lc rgb "#FF6347", "CorrelationINF.txt"  u (((abs($1)))*(1/(1))):(((abs($2)))*(1/(1))) t 'D {/Symbol \256 } {/Symbol \245 }' with points  pointtype 6 lw 3  ps 5.0 lc rgb "#4682B4", 


#M(x) t 'C(r)=r^{-2.9}' with line dt 4 lw 4 linetype 4 lc rgb "blue"



