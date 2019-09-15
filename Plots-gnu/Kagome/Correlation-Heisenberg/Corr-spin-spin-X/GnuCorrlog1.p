set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"
set key t r
set key box 6
set key font ",38"
set key spacing 5
set key width -1


M4(x)=c4*exp(-x*d4)
fit  [20:60] M4(x) "CorrelationHD6-two.txt" u ($1):(abs($2))    via c4,d4 

M5(x)=c5*exp(-x*d5)
fit [10:60]  M5(x) "CorrelationH-two.txt" u ($1):(abs($2))    via c5,d5 

#M6(x)=c6*exp(-x*d6)
#fit  [15:40] M6(x) "CorrelationH7.txt" u ($1):(abs($2))    via c6,d6 

#M7(x)=c7*exp(-x*d7)
#fit  [12:50] M7(x) "CorrelationH8.txt" u ($1):(abs($2))   via c7,d7 


#M8(x)=c8*exp(-x*d8)
#fit  [12:50] M8(x) "CorrelationH9.txt" u ($1):(abs($2))    via c8,d8 


set tics scale 4
set format y ""
#set format x ""
#set logscale x 10
set logscale y 10
set ytics 1e-13, 10, 1
#set xtics 1e-14, 10, 1
#set mxtics 10
set ytics add ("10^{-1}" 0.1,"10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001,"10^{-5}" 0.00001, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001,"10^{-9}" 0.000000001,"10^{-10}" 0.0000000001,"10^{-11}" 0.00000000001,"10^{-11}" 0.00000000001,"10^{-12}" 0.000000000001,"10^{-13}" 0.0000000000001)
set xtics add ("10^{0}" 1.0, "10^{-1}" 0.1, "10" 10,"10^{-2}" 0.01, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)
set mytics 14
set xtics nomirror
show logscale

set output "Corrlog1.eps"

p  [2:60] [0.000000000001:0.01] "CorrelationHD6-two.txt"  u ($1):(abs($2)) t 'two-layer: D=6'  with points  pointtype 9 lw 3  ps 8 lc rgb "#DC143C","CorrelationHD6-one.txt"  u ($1):(abs($2)) t 'one-layer: D=6'   with points  pointtype 8 lw 3  ps 8 lc rgb "#008080","CorrelationH-two.txt" u ($1):(abs($2)) t 'two-layer: D=5'  with points  pointtype 7 lw 3  ps 8 lc rgb "#008080","CorrelationH-one.txt" u ($1):(abs($2)) t 'one-layer: D=5'  with points  pointtype 6 lw 3  ps 8 lc rgb "#FF1493",M4(x) notitle with line dt 5 lw 4 linetype 4 lc rgb "#191970",M5(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970"