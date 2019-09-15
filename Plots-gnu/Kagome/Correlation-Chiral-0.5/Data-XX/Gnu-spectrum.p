set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "{/Symbol f}"   font ",40" textcolor rgb "black"
#set ylabel "-{/Symbol l}"  font ",40"  textcolor rgb "black"

#set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);


#set key font ",35"
#set key spacing 8
#set key samplen 6

#set key box vertical width 0.5 height 1 samplen 6 b l
#set key box b l  

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "spectrum.eps"
p  [-3.0:3.0] [0:5]  "Spectrum.txt"  u 3:4  notitle   with points  pointtype 6 lw 2 ps 6.0 lc rgb "#4169E1"




