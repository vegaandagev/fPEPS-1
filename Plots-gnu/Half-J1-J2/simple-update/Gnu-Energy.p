set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",40" textcolor rgb "black"
set ylabel '{/Symbol D}E'  font ",40"  textcolor rgb "black"
#set format y "10^{%L}"
#set ytics(0.01,0.001,0.0001,0,00001)
set xtics(2,6,10,14)
set logscale y
#set ytics font ", 30"
set tics scale 4
#set mxtics 4
#set mytics 4
set format y ""
set ytics 1e-6, 10, 1
set ytics add ("1" 1, "10^{-2}" 0.01, "10^{-3}" 0.001, "10^-4" 0.0001)
set mytics 10
set key box t r
set key font ",40"
set key spacing 8  
set key font ",40"
set key spacing 8
set output "HeisenbergSU.eps"
p [0:17] [0.0007:0.015] "Nosymmetry.txt" u ($1):(($2)*(1/(0.66944))) t  "No symmetry"  with linespoints lw 3 lt rgb "#A52A2A" pointtype 2  ps 7.0, "Usymmetry.txt" u ($1):(($2)*(1/(0.66944))) t "U(1)" with linespoints lt rgb "#3CB371" lw 3 pointtype 4  ps 7.0, "previous.txt"  t "Bauer et al." with linespoints lw 3  pointtype 6  ps 7.0
