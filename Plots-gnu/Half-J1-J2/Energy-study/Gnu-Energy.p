set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "iteration"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"
#set format y "10^{%L}"
set ytics(-0.494,-0.4945,-0.4935)
#set logscale y
set tics scale 4
set key box t r
set key font ",35"
set key spacing 6  
set key font ",40"
set key spacing 6
set output "HeisenbergFU.eps"
p [0:35]  "CG.txt"  title "full CG" with linespoints lw 3 pointtype 2  ps 6 lt rgb "#A52A2A", "SVD.txt" t "New Scheme" with linespoints lw 3 pointtype 14  ps 6.0  lt rgb "#2E8B57", "SVD-mpo.txt"  title "Corboz et al" with linespoints lw 3 pointtype 8  ps 6 lt rgb "#DC143C"








