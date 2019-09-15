set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "N_{v}"   font ",40" textcolor rgb "black"
#set ylabel "S"  font ",40"  textcolor rgb "black"

#set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);

G(x)=b+a*x
fit   G(x) "S.txt" u ($1):($2)    via a,b

set key font ",35"
set key spacing 8
set key samplen 6

set key box vertical width 0.5 height 1 samplen 6 b l
set key box t l  

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "VNE.eps"
p [0:9]   "S.txt" u ($1):(abs($2)) notitle  with points  pointtype 2 lw 2  ps 3.0 lc rgb "#8B008B",x<9.0 ? G(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4





