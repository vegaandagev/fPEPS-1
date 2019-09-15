set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"

M(x)=d1*x**c1
fit [4:12]  M(x) "CorrelationHX.txt" u ($1):(abs($2))   via c1,d1 

#M(x)=(c1*x)+d1

#fit  [log(1)/log(10):log(7)/log(10)] M(x) "CorrelationHF.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    via #c1,d1 

set key b l 
set key box 6
set key font ",40"
set tics scale 4
set format y ""
set format x ""
set logscale x 10
set logscale y 10
set ytics 1e-10, 10, 1
#set xtics 1e-10, 10, 1
set mytics 10
set mxtics 10
set ytics add ("10^{-1}" 0.1,"10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001,"10^{-5}" 0.00001, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)
set xtics add ("10^{0}" 1.0, "10^{-1}" 0.1, "10^{1}" 10,"10^{-2}" 0.01, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)

show logscale

set key spacing 6
set output "Corrlog.eps"

p  [1.0:45]   "Corr.txt"  u ($1):(abs($2)) t 'spin:D=10'  with points  pointtype 2 lw 3  ps 8 lc rgb "#FF1493","CorrelationHDimer.txt"  u ($1):(abs($2)) t 'Dimer:D=10'  with points  pointtype 4 lw 3  ps 8 lc rgb "#9400D3","CorrelationHChiral.txt" u ($1):(abs($2)) t 'chiral:D=10'  with points  pointtype 6 lw 3  ps 8 lc rgb "#008080"


