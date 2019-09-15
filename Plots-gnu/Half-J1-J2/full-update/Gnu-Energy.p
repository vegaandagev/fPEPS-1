set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",40" textcolor rgb "black"
#set ylabel '{/Symbol D}E'  font ",40"  textcolor rgb "black"
set logscale y
#set ytics font ", 30"
set tics scale 4
#set mxtics 4
#set mytics 4
#set format y ""
#set ytics 1e-6, 10, 1
set ytics add ("1" 1, "10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001)
#set mytics 10

set key box t r
set key font ",40"
set key spacing 8
set output "HeisenbergFU.eps"
p [1.4:7.9] [0.00008:0.015] "Usymmetry.txt"  t "U(1)" with linespoints lt rgb "#3CB371" lw 3  pointtype 4  ps 7.0, "Nosymmetry.txt"  t "No-symmetry" with linespoints lt rgb "#A52A2A" lw 3 pointtype 2  ps 7


