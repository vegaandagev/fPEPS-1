set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "Time[second]"   font ",35" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"
#set ytics nomirror (-0.6682, -0.6690, -0.6692, -0.6694)
set ytics 0.00025
set key box t r  
set key font ",35"
set key spacing 1.5
set key width -6

set output "EnergyComparing.eps"

p [6:13] [-0.66943: -0.6682] "Ef.txt" u (log($1)):($2) title "two-layer: (D, {/Symbol c})=(5,25)"  with points  pointtype 2 lw 2  ps 8.0 lc rgb "#8B008B", "Eone.txt" u (log($1)):($2) title "one-layer: (D, {/Symbol c})=(5,40)"  with points  pointtype 4 lw 2  ps 8.0 lc rgb "#a40000", -0.6694 t "Monte Carlo" with line lw 4 linetype 7 lc rgb "#2F4F4F"
