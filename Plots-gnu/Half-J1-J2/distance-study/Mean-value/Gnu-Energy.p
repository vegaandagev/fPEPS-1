set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "u"   font ",45" textcolor rgb "black"
set ylabel ' '  font ",45"  textcolor rgb "black"
#set format y "10^{%L}"
#set ytics(1.5,1.25,1.0,0.75,0.5)
set tics scale 4
set logscale y
set format y ""
set ytics 1e-10, 10, 1
set ytics add ("10^{-2}" 0.01, "10^{-3}" 0.001, "10^{-4}" 0.0001,"10^{-5}" 0.00001, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)
set mytics 10
set key box t r
set key font ",40"
set key spacing 8
set output "HeisenbergFU.eps"
p  [1:30]   [0.000000001:0.01] "SVD.txt" u (($1)):(($3)/40174.436671) t "New scheme" with linespoints lw 3 pointtype 12  ps 6.0 lt rgb "#FF1493", "SVD-mpo.txt" u ((($1)-10)):(($3)/40174.436671) t "Corboz et al" with linespoints lw 3 pointtype 8  ps 6.0 lt rgb "#DC143C"

#, "SVDmix.txt" u ((($1))):((($3))/40174.436671) t "steps-(a-iv)" with linespoints lw 3 pointtype 14  ps 6.0  lt rgb "#2E8B57"
#[0.62:1.5]

#40174.436671