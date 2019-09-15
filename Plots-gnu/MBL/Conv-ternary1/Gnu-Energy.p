set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "iteration"   font ",40" textcolor rgb "black"
set ylabel ' {/Symbol s} / 2^{N}'  font ",40"  textcolor rgb "black"
#set format y "10^{%L}"
#set ytics(1.5,1.25,1.0,0.75,0.5)
set tics scale 4
set logscale y 10
set logscale x 10
set format y ""
#set format x ""
set ytics 1e-7, 10, 10
set xtics 1e+7, -10, 1
set ytics add ("10" 10,"1" 1, "10^{-1}" 0.1,"10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001,"10^{-5}" 0.00001,"10^{-6}" 0.000001)
set xtics add ("10" 10,"1" 1,"10^{2}" 100,"10^{3}" 1000,"10^{5}" 100000,"10^{4}" 10000, "10^{6}" 1000000, "10^{8}" 100000000)
set mytics 10
set mxtics 10
#set ytics nomirror
set xtics nomirror
#set mxtics 10
set key box t r
#set key outside
set key font ",35"
set key spacing 1
set key width -3
#set key at 1000000.5,12.9
set output "Cost.eps"
p [100:90000] [0.07:8]  "varianceSVD_0.txt"  u ($1):($2) t "Linearizing" with points lw 3 lt rgb "#FF0000" pointtype 2  ps 6, "variancePol_0.txt" u ($1):($2) t "CG-Poly" with points lw 3 pointtype 4  ps 6  lt rgb "#FF1493","varianceArm_0.txt" u ($1):($2) t "CG-Armji" with points lw 3 pointtype 8  ps 6 lt rgb "#6A5ACD" 


#40174.436671
