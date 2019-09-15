#set terminal postscript eps enhanced color font 'Helvetica,14' lw 3
set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"
set key b l 
set key box 6
set key font ",40"

#set key width 6

set key spacing 6


M4(x)=c4* exp(-x*d4)
fit  [10:30] M4(x) "CorrelationH4.txt" u ($1):(abs($2))    via c4,d4 

M5(x)=c5* exp(-x*d5)
fit  [10:30] M5(x) "CorrelationH5.txt" u ($1):(abs($2))    via c5,d5 

M6(x)=c6*exp(-x*d6)
fit  [10:30] M6(x) "CorrelationH6.txt" u ($1):(abs($2))    via c6,d6 

M7(x)=c7*exp(-x*d7)
fit  [10:30] M7(x) "CorrelationH7.txt" u ($1):(abs($2))   via c7,d7 


M8(x)=c8*exp(-x*d8)
fit  [10:30] M8(x) "CorrelationH8.txt" u ($1):(abs($2))    via c8,d8 



set tics scale 4
set format y ""
#set format x ""
#set logscale x 10
set logscale y 10
set ytics 1e-14, 10, 1
#set xtics 1e-10, 10, 1
#set mxtics 10
set ytics add ("10^{-2}" 0.01,"10^{-4}" 0.0001, "10^{-6}" 0.000001, "10^{-8}" 0.00000001,"10^{-10}" 0.0000000001, "10^{-12}" 0.000000000001)
set xtics add ("10^{0}" 1.0, "10^{-1}" 0.1, "10" 10,"10^{-2}" 0.01, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)
set mytics 14
set xtics nomirror
show logscale


set output "Corrlog1.eps"


#"CorrelationH4.txt"  u ($1):(abs($2)) t 'D=4'  with points  pointtype 6 lw 3  ps 5.0 lc rgb "#48D1CC",

p  [5:40]  "CorrelationH4.txt"  u ($1):(abs($2)) t 'D=4'  with points  pointtype 6 lw 3  ps 5.0 lc rgb "#DC143C", "CorrelationH5.txt"  u ($1):(abs($2)) t 'D=5'  with points  pointtype 5 lw 3  ps 5.0 lc rgb "#FF69B4", "CorrelationH6.txt"  u ($1):(abs($2)) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationH7.txt"  u ($1):(abs($2)) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED", "CorrelationH8.txt"  u ($1):(abs($2)) t 'D=8' with points pointtype 12 lw 3  ps 5.0 lc rgb "#FF6347", M4(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",M5(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", M6(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",  M7(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",  M8(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970"

#M4(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970"

