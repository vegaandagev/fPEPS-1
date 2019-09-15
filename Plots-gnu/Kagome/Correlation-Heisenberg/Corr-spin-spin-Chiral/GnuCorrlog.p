set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C_{t}^{s}'  font ",40"  textcolor rgb "black"

M(x)=d1*x**c1
fit [4:12]  M(x) "CorrelationH.txt" u ($1):(abs($2))   via c1,d1 

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

p  [1.0:25]   "CorrelationH.txt"  u ($1):(abs($2)) t 'D=10'  with points  pointtype 2 lw 3  ps 8 lc rgb "#4682B4"

#,"CorrelationV6.txt"  u ($1):(abs($2)) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371"#, "CorrelationV7.txt"  u ($1):(abs($2)) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED", "CorrelationV8.txt"  u ($1):(abs($2)) t 'D=8' with points  pointtype 12 lw 3  ps 5.0 lc rgb "#FF6347", "CorrelationV9.txt"  u ($1):(abs($2))  t 'D=9' with points  pointtype 14 lw 3  ps 5.0 lc rgb "#9400D3", "CorrelationV12.txt"  u ($1):(abs($2))  t 'D=12' with points  pointtype 16 lw 3  ps 5.0 lc rgb "#8f5902" ,  M(x) t 'r^{-1.6}' with line dt 4 lw 4 linetype 4 lc rgb "#FF1493"

#, "CorrelationH5.txt"  u ($1):(abs($2)) t 'D=4'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#8B0000","CorrelationHH4.txt"  u ($1):(abs($2)) t 'D=4'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#9400D3"

#,
#,((-1.75)*x)+d1 t 'r^{-1.7}' with line dt 4 lw 4 linetype 4 lc rgb "red"

#,((-1.55)*x)+d1 t 'r^{-1.5}' with line dt 4 lw 4 linetype 4 lc rgb "blue”




