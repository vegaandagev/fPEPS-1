set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "D"   font ",40" textcolor rgb "black"
set ylabel "N_{CTM}"  font ",40"  textcolor rgb "black"

#set xtics ("1/1000" 0.001,"1/250" 0.004, "1/100" 0.01)
#set xtics (0.3,0.2,0.1);
#set ytics 0.0001



set key box t l  
set key font ",35"
set key spacing 1.5

set output "Niteration.eps"

p [1.9:8.1] [0:62] "CTM-sppedup.txt" u ($1):($2) title "two-layer"  with points  pointtype 2 lw 2  ps 8.0 lc rgb "#8B008B", "CTM-sppedup.txt" u ($1):($3) title "one-layer"  with points  pointtype 4 lw 2  ps 8.0 lc rgb "#a40000"






