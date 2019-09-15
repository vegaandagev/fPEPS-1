set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "W"   font ",40" textcolor rgb "black"
set ylabel ' {/Symbol s} / 2^{N}'  font ",40"  textcolor rgb "black"
#set format y "10^{%L}"
#set ytics(1.5,1.25,1.0,0.75,0.5)
set tics scale 4
set logscale y 10
#set logscale x 10
set format y ""
#set format x ""
set ytics 1e-7, 10, 30
#set xtics 1e+7, -10, 1
set ytics add ("10" 10,"1" 1, "10^{-1}" 0.1,"10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001,"10^{-5}" 0.00001,"10^{-6}" 0.000001)
set xtics add (2,4,6,8,10,12,14)
set mytics 10
#set mxtics 10
#set ytics nomirror
set xtics nomirror
#set mxtics 10
set key box b l
#set key outside
set key font ",36"
set key spacing 6
set key width 1
#set key at 1000000.5,12.9
set output "Cost.eps"
p [1:10] [0.005:1] "varianceAll.txt"  u ($1):($2) t "l=6" with points lw 3 lt rgb "#FF1493" pointtype 2  ps 8,"varianceAll.txt"  u ($1):($3) t "l=4" with points lw 3 lt rgb "#9400D3" pointtype 4  ps 8, "varianceAll.txt"  u ($1):($4) t "l=2" with points lw 3 lt rgb "#FFA500" pointtype 6  ps 8, "varianceAll.txt"  u ($1):($6) t "Binary" with points lw 3 lt rgb "#A0522D" pointtype 12  ps 8



#"varianceAll.txt"  u ($1):($5) t "Wahl:l=6" with points lw 3 lt rgb "#DC143C" pointtype 8  ps 8
#40174.436671