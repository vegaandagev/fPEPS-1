set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol x}"   font ",40" textcolor rgb "black"

set key b l 
set key box 6
set key font ",40"
set key spacing 6
set output "Length.eps"


set xtics (0,5, 6,7,8)
#set ytics nomirror ( -3, -5, -6, -14);





#fit  M(x) "L.txt" u ((log(abs($1)))*(1/1)):((log(abs($2)))*(1/1))    via c1,d1 

M(x)=(c1**x)+d1

fit  M(x) "L.txt"    via c1,d1 



#p   "L.txt"  u ((log(abs($1)))*(1/1)):((log(abs($2)))*(1/1)) t "Corr"   with points  pointtype 2 lw 3  ps 5.0 lc rgb "#4B0082", M(x) t 'C(r)=r^{-2.9}' with line dt 4 lw 4 linetype 4 lc rgb "blue"
p  [0:7] "L.txt"  t "Corr"   with points  pointtype 2 lw 4  ps 6.0 lc rgb "#FF1493", M(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2"


