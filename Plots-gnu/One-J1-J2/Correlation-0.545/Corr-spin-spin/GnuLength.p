set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",48" textcolor rgb "black"
set ylabel "{/Symbol x}"   font ",48" textcolor rgb "black"

set key b r 
set key box 4
set key font ",38"
set key spacing 8
set output "Length.eps"
set tics scale 3
set tics scale 3


set xtics nomirror (0,4,5, 6,7,8)
set ytics nomirror ( 0,1, 2, 3, 4);


M(x)=(x**(c1))*d1

fit [1:9] M(x) "L1.txt"    via c1,d1 

M2(x)=(x**(c2))*d2

fit  M2(x) "L1.txt" u ($1):($3)   via c2,d2 

#M3(x)=(x**(c3))*d3

#fit  M3(x) "L1.txt" u ($1):($4)   via c3,d3 


#p   "L.txt"  u ((log(abs($1)))*(1/log(2))):((log(abs($2)))*(1/log(2))) t "Corr"   with points  pointtype 2 lw 3  ps 5.0 lc rgb "#4B0082", M(x) t 'C(r)=r^{-2.9}' with line dt 4 lw 4 linetype 4 lc rgb "blue"


p  [3.5:8.5] [0.9:2] "L1.txt" u  ($1):($2) t "spin-spin"   with  points  pointtype 2 lw 3  ps 8.0 lc rgb "#FF1493","L1.txt" u  ($1):($3) t "spin-spin"   with  points  pointtype 2 lw 3  ps 8.0 lc rgb "#FF1493", M(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2", M2(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#483D8B"


#,M3(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#6A5ACD"
#, "L1.txt" u  ($1):($3) t "dimer-dimer"  with  points  pointtype 6 lw 3  ps 8.0 lc rgb "#EE82EE", "L1.txt" u  ($1):($4) t "transfer-matrix"  with  points  pointtype 4 lw 3  ps 8.0 lc rgb "#008080",

