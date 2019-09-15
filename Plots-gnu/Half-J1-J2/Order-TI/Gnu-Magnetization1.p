set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 55”
set ytics font "Times-Roman, 55”

set xlabel "1/D"   font ",55" textcolor rgb "black"
set ylabel "{/Symbol D}T_{x}"  font ",55"  textcolor rgb "black"
#set log y 
set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
set ytics (0.4,0.2,0.3,0.1,0.05,0.01)




#set key box 2
#set key c r 
#set key font ",30"
#set key spacing 6

#set key at 0.1,0.7
set output "TI1.eps"
p  [0.008:0.35] [0.0:0.10]  'J20.6.txt' u ($1):($3)  notitle 'J2=0.60'  with points  pointtype 3 lw 3  ps 6.0 lc rgb "#FF1493",'J20.55.txt' u ($1):($3)  notitle 'J2=0.55'   with points  pointtype 2 lw 3  ps 6.0 lc rgb "#8B008B", 'J20.54.txt' u ($1):($3)  notitle 'J2=0.54'  with points  pointtype 4 lw 3  ps 6.0 lc rgb "#8A2BE2"
