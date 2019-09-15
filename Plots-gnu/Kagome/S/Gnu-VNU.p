set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "Length"   font ",40" textcolor rgb "black"
#set ylabel "S"  font ",40"  textcolor rgb "black"

#set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);

#G(x)=b+a*x
#fit   G(x) "S.txt" u ($1):($2)    via a,b

M(x)=((x**a1))*b1
fit [4:40]  M(x) "SpectrumC4.txt" u ($1):($2)    via a1,b1

G(x)=((x**a))*b
fit [10:45]  G(x) "SpectrumC6.txt" u ($1):($2)    via a,b

#P(x)=((x**a2))*b2
#fit [20:40]  P(x) "SpectrumC8.txt" u ($1):($2)    via a2,b2

set key font ",35"
set key spacing 8
set key samplen 6

#set key box t l  

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "SvsC.eps"
#p  [20:40] "SpectrumC8.txt" u ($1):($2) title "w=8"  with points  pointtype 6 lw 2  ps 8.0 lc rgb "#C71585", x<70 ? P(x) : 1/0 notitle lw 4 lc rgb "#a40000" dt 4

p  [1:100] "SpectrumC6.txt" u ($1):($2) title "w=6"  with points  pointtype 6 lw 2  ps 8.0 lc rgb "#C71585", x<70 ? G(x) : 1/0 notitle lw 4 lc rgb "#a40000" dt 4


p  [1:100] "SpectrumC4.txt" u ($1):($2) title "w=4"  with points  pointtype 2 lw 2  ps 8.0 lc rgb "#8B008B", x<70 ? M(x) : 1/0 notitle lw 4 lc rgb "#a40000" dt 4





