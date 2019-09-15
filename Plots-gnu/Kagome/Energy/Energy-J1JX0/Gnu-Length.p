set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol x}"  font ",40"  textcolor rgb "black"

set xtics nomirror ("3" 3, "5" 5, "8" 8,"11" 11)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key b r 
set key box 4
set key font ",38"
set key spacing 8
set output "Length.eps"
set tics scale 3
set tics scale 3
set key width 4

M(x)=(x**(c1))*d1

fit  M(x) "EKagome.txt"  u (1/$1):($4)  via c1,d1 

G(x)=(x**(c2))*d2

fit  G(x) "Echiral5.txt"  u (1/$1):($4)  via c2,d2



set output "Length.eps"
p [0.0:11.0]  [0:2]  "EKagome.txt"  u (1/$1):($4)  t "J_{{/Symbol c}}=0"  with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B","Echiral5.txt"  u (1/$1):($4)  t "J_{{/Symbol c}}=0.5"  with points  pointtype 4 lw 3  ps 8.0 lc rgb "#ce5c00", M(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#f57900"#, G(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2"

#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,

