set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "iteration"   font ",45" textcolor rgb "black"
set ylabel ' {/Symbol s} / 2^{N}'  font ",45"  textcolor rgb "black"
#set format y "10^{%L}"
#set ytics(1.5,1.25,1.0,0.75,0.5)
set tics scale 4
set logscale y 10
set logscale x 10
set format y ""
#set format x ""
set ytics 1e-6, 10, 10
set xtics 1e+6, -10, 1
set ytics add ("10" 10,"1" 1, "10^{-1}" 0.1,"10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001,"10^{-5}" 0.00001)
set xtics add ("10" 10,"1" 1,"10^{2}" 100,"10^{3}" 1000,"10^{5}" 100000,"10^{4}" 10000, "10^{6}" 1000000, "10^{8}" 100000000)
set mytics 10
set mxtics 10
set ytics nomirror
set xtics nomirror
#set mxtics 10
set key box t r
set key font ",40"
set key spacing 8
set output "Cost.eps"
p [20:200000] [0.2:20]  "variance1.txt"  u ($1):($2) t "Linearizing" with points lw 3 lt rgb "#A52A2A" pointtype 2  ps 6.0, "variance2.txt" u ($1):($2) t "Poly CG" with points lw 3 pointtype 4  ps 6.0  lt rgb "#2E8B57" ,"variance3.txt" u ($1):($2) t "Armji CG" with points lw 3 pointtype 6  ps 6.0 lt rgb "#FF1493" , "variance4.txt" u ($1):($2) t "Armji SD" with points lw 3 pointtype 8  ps 6.0 lt rgb "#DC143C","variance5.txt" u ($1):($2) t "Poly SD" with points lw 3 pointtype 10  ps 6.0 lt rgb "#DC143C"



#40174.436671