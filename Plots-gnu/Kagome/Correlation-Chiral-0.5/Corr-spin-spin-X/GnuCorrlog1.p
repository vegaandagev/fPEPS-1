set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C(r)'  font ",40"  textcolor rgb "black"
set key t r 
set key box 6
set key font ",40"
set key spacing 6


M4(x)=c4*exp(-x*d4)
fit  [15:60] M4(x) "CorrelationHZ.txt" u ($1):(abs($2))    via c4,d4 

M5(x)=c5*exp(-x*d5)
fit [10:20]  M5(x) "CorrelationHChiral.txt" u ($1):(abs($2))    via c5,d5 


M6(x)=c6*exp(-x*d6)
fit [5:20]  M6(x) "CorrelationHDimerE.txt" u ($1):(abs($2))    via c6,d6 


set tics scale 4
set format y ""
#set format x ""
#set logscale x 10
set logscale y 10
set ytics 1e-14, 10, 1
#set xtics 1e-10, 10, 1
#set mxtics 10
set ytics add ("10^{-2}" 0.01,"10^{-4}" 0.0001, "10^{-6}" 0.000001, "10^{-8}" 0.00000001,"10^{-10}" 0.0000000001, "10^{-12}" 0.000000000001)
set xtics add ("1" 1.0, "10^{-1}" 0.1, "10" 10,"10^{-2}" 0.01, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)
set mytics 14
set xtics nomirror
show logscale


set output "Corrlog1.eps"


p  [1:30] [0.00000000001:0.001] "CorrelationHZ.txt" u ($1):(abs($2)) t 'spin: D=10' with points  pointtype 15 lw 3 ps 8 lc rgb "#FF1493", "CorrelationHDimerE.txt" u ($1):(abs($2)) t 'dimer: D=10'  with points  pointtype 17 lw 3  ps 8 lc rgb "#204a87",M4(x) notitle with line dt 5 lw 4 linetype 4 lc rgb "#191970",M6(x) notitle with line dt 5 lw 4 linetype 4 lc rgb "#191970"


#, M5(x) notitle with line dt 5 lw 4 linetype 4 lc rgb "#ce5c00"

#"CorrelationHDimer.txt" u ($1):(abs($2)) t 'dimer: D=10'  with points  pointtype 4 lw 3  ps 8 lc rgb "#008080"

#"CorrelationHChiral.txt"u ($1):(abs($2)) t 'chiral: D=10'  with points  pointtype 6 lw 3  ps 8 lc rgb "#008080"