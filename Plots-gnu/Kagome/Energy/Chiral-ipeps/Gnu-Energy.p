set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "J_{c}"   font ",40" textcolor rgb "black"
set ylabel "H_c"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key box b r  
set key font ",40"
set key spacing 6
set key width 2



G1(x)=0

set output "chiral.eps"
set multiplot
p [0:0.30] [0:0.25] "M7.txt" u ($1):($2)  t "D=7" with points   pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C", "M9.txt" u ($1):($2)  t "D=9" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#C71585", "M11.txt" u ($1):($2)  t "D=11"  with points  pointtype 17 lw 3  ps 8.0 lc rgb "#8f5902","chiral.txt" u ($1):(abs($2))  t "DMRG" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#FF69B4"

#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,


set origin .2, .60
set size 0.55,0.320
clear
unset key
unset xlabel
set ylabel "{/Symbol \266}H_{c}/{/Symbol \266}J_{c}"  font ",30"  textcolor rgb "black"
set ytics font "Times-Roman, 30”
set xtics font "Times-Roman, 30”
set tics scale 2
set xtics nomirror

set ytics (0.4,0.2,0.3,0.5,0.6,0.7)


p [.0:.26]    "Der7.txt" u ($1):(abs($2))  t "D=7" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C","Der9.txt" u ($1):(abs($2))  t "D=9" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#C71585","Der11.txt" u ($1):(abs($2))  t "D=11" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#8f5902"#,"DMRG.txt" u ($1):(abs($2))  t "DMRG" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#DA70D6"

unset multiplot






