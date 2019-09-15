set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
#set x2tics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/{/Symbol c}"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol x}"  font ",40"  textcolor rgb "black"
#set x2label "1/{/Symbol c}_{one-layer}" font ",40" textcolor rgb "black"

set tics scale 3
set ytics 0.5
set xtics 0.01
#set x2tics 0.02

G(x)=b+a*x
fit   G(x) "L4.txt" u (1/$1):($2)    via a,b

G1(x)=b1+a1*x
fit   G1(x) "L5.txt" u (1/$1):($2)    via a1, b1

G2(x)=b2+a2*x
fit   G2(x) "L6.txt" u (1/$1):($2)    via a2, b2

#set key samplen 4
#set key box vertical width 0.5 height 1 samplen 6 b l
set key box t r  
set key font ",36"
set key spacing 1

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "Length.eps"

#set multiplot
p [0:0.03] [1:4.1] "L4.txt" u (1/$1):($2) title "D=4: one-layer" axes x1y1 with points  pointtype 5 lw 2  ps 8.0 lc rgb "#a40000","L5.txt" u (1/$1):($2) title "D=5: one-layer"  with points  pointtype 7 lw 2  ps 8.0 lc rgb "#1e90ff","L6.txt" u (1/$1):($2) title "D=6: one-layer"  with points  pointtype 9 lw 2  ps 8.0 lc rgb "#8B008B",x<0.015 ? G(x) : 1/0 notitle lw 4 lc rgb "#8B008B" dt 4,x<0.015 ? G1(x) : 1/0 notitle lw 4 lc rgb "#a40000" dt 4, x<0.015 ? G2(x) : 1/0 notitle lw 4 lc rgb "#5c3566" dt 4,"L4F.txt" u (1/$1):($2) title "D=4: two-layer" with points  pointtype 4 lw 2  ps 8.0 lc rgb "#2E8B57","L5F.txt" u (1/$1):($2) title "D=5: two-layer"  with points  pointtype 6 lw 2  ps 8.0 lc rgb "#DC143C","L6F.txt" u (1/$1):($2) title "D=6: two-layer" with points  pointtype 8 lw 2  ps 8.0 lc rgb "#FF6347"


#set key box t l  

#plot [0:0.01] [1:4.1] "L4F.txt" u (1/$1):($2) title "D=4" axes x2y1 with points  pointtype 4 lw 2  ps 8.0 lc rgb "#2E8B57","L5F.txt" u (1/$1):($2) title "D=5" axes x2y1 with points  pointtype 6 lw 2  ps 8.0 lc rgb "#DC143C","L6F.txt" u (1/$1):($2) title "D=6" axes x2y1 with points  pointtype 8 lw 2  ps 8.0 lc rgb "#FF6347"
