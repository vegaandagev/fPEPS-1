set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "l"   font ",35" textcolor rgb "black"
set ylabel 'Time[s]'  font ",35"  textcolor rgb "black"
set tics scale 4
set logscale y 10
set format y ""
set ytics 1e+7, -10, 1
set ytics add ("10" 10,"1" 1, "10^{1}" 10,"10^{2}" 100.0, "10^{3}" 1000, "10^{4}" 10000,"10^{5}" 100000,"10^{6}" 1000000,"10^{7}" 10000000)
set xtics add (2,4,6)
set mytics 10
set key box b r
#set key outside
set key font ",35"
set key spacing 2
set key width 2
set output "Cost.eps"
set xtics nomirror

set multiplot

p [0.5:6] [4:2000000] "time.txt"  u ($1):($2)*(10000) t "CG-Poly" with linespoints lw 2.5 lt rgb "#EE82EE" pointtype 4  ps 8, "time.txt"  u ($1):($3)*(10000) t "CG-Armji" with linespoints lw 2.5 lt rgb "#6A5ACD" pointtype 8  ps 8

set origin .15, .56
set size 0.45,0.4

clear
unset key
unset xlabel
unset ylabel
unset logscale y
unset format y 
unset ytics 
unset mytics 

set xtics add (2,4,6)

set ylabel "speed-up factor"  font ",35"  textcolor rgb "black"

set xtics font "Times-Roman, 25”
set ytics font "Times-Roman, 25”

set tics scale 1.5

set xtics (2,4,6)
set ytics add (4,6,8,10,12)

set xtics nomirror
#set ytics nomirror


p "time.txt" u ($1):($4)  notitle  with linespoints  pointtype 2 lw 2  ps 6.0 lc rgb "#FF1493"

unset multiplot

