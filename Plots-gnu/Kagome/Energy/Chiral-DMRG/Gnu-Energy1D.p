set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "J_{ch}"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol \266}<H_{c}>/{/Symbol \266}J_{ch}"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key box b l  
set key font ",38"
set key spacing 6
set key width 2
G1(x)=0

set output "E-D.eps"
set multiplot

p [-.01:.28] [0.1:1.2]    "Der7.txt" u ($1):(abs($2))  t "D=7" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C","Der9.txt" u ($1):(abs($2))  t "D=9" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#C71585","Der11.txt" u ($1):(abs($2))  t "D=11" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#8f5902","DMRG.txt" u ($1):(abs($2))  t "DMRG" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#DA70D6"
 
 
 
set origin .36, .540
set size 0.58,0.40
clear
unset key
unset xlabel
set ylabel "{/Symbol \266}^{2}<H_{c}>/{/Symbol \266}{J_{ch}}^{2}"  font ",30"  textcolor rgb "black"
set ytics font "Times-Roman, 30”
set xtics font "Times-Roman, 30”
set tics scale 2
set xtics nomirror

set xtics (0.1,0.2,0.3)
set ytics (0.4,0.8,1.2,1.6,2)


p [.05:.35]    "Der7.txt" u ($1):(abs($3))  notitle with points  pointtype 9 lw 3  ps 3.0 lc rgb "#DC143C","Der9.txt" u ($1):(abs($3))  notitle with points  pointtype 13 lw 3  ps 3.0 lc rgb "#C71585","Der11.txt" u ($1):(abs($3))  notitle with points  pointtype 17 lw 3  ps 3.0 lc rgb "#8f5902"

unset multiplot
 
