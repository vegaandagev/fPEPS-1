set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",48" textcolor rgb "black"
set ylabel "{/Symbol x}"   font ",48" textcolor rgb "black"

set key t l 
set key box 4
set key font ",40"
set key spacing 8
set output "Length.eps"
set tics scale 4


set xtics nomirror (4,5, 6,7,8)
set ytics nomirror ( 3, 2, 1, 0);





#fit  M(x) "L.txt" u ((log(abs($1)))*(1/log(2))):((log(abs($2)))*(1/log(2)))    via c1,d1 

#M(x)=(x**(c1))*d1

#fit  M(x) "L.txt"    via c1,d1 

#M2(x)=(x**(c2))*d2

#fit  M2(x) "L.txt" u ($1):($3)   via c2,d2 

#M3(x)=(x**(c3))*d3

#fit  M3(x) "L.txt" u ($1):($4)   via c3,d3 

#M4(x)=(x**(c4))*d4

#fit  M4(x) "L.txt" u ($1):($4)   via c4,d4 

M5(x)=(x**(c5))*d5

fit  M5(x) "L.txt" u ($1):($6)   via c5,d5 

M6(x)=(x**(c6))*d6

fit  M6(x) "L.txt" u ($1):($7)   via c6,d6 


p  [3.0:9.1] [1:3.2] "L.txt" u  ($1):($6) t "transfer-matrix {/Symbol x}_{y}"   with  points  pointtype 2 lw 3  ps 8.0 lc rgb "#FF1493", "L.txt" u  ($1):($7) t "transfer-matrix, {/Symbol x}_{x}"  with  points  pointtype 6 lw 3  ps 8.0 lc rgb "#3CB371", M5(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2", M6(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#483D8B"

#, M(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2", M2(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#483D8B", M3(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2", M4(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#483D8B"

#"L.txt" u  ($1):($2) t "spin-spin, {/Symbol x}_{y}"   with  points  pointtype 2 lw 3  ps 8.0 lc rgb "#FF1493", "L.txt" u  ($1):($3) t "spin-spin, {/Symbol x}_{x}"  with  points  pointtype 3 lw 3  ps 8.0 lc rgb "#3CB371","L.txt" u  ($1):($4) t "dimer-dimer, {/Symbol x}_{y}"   with  points  pointtype 6 lw 3  ps 8.0 lc rgb "#FF1493", "L.txt" u  ($1):($5) t "dimer-dimer, {/Symbol x}_{x}"  with  points  pointtype 7 lw 3  ps 8.0 lc rgb "#3CB371",

