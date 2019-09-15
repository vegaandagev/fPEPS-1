set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/{/Symbol c}"   font ",35" textcolor rgb "black"
set ylabel "{/Symbol x}"  font ",40"  textcolor rgb "black"

#set xtics ("1/1000" 0.001,"1/250" 0.004, "1/100" 0.01)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);

G(x)=b+a*x
fit [0:0.0025]  G(x) "L6.txt" u (1/$1):($2)    via a,b

G1(x)=b1+a1*x
fit  [0:0.0055] G1(x) "L7.txt" u (1/$1):($2)    via a1, b1

G2(x)=b2+a2*x
fit  [0:0.0055] G2(x) "L8.txt" u (1/$1):($2)    via a2, b2

G3(x)=b3+a3*x
fit   G3(x) "L9.txt" u (1/$1):($2)    via a3, b3

G4(x)=b4+a4*x
fit [0:0.0055] G4(x) "L13.txt" u (1/$1):($2+0.09)    via a4, b4

G5(x)=b5+a5*x
fit  G5(x) "L12.txt" u (1/$1):($2+0.09)    via a5, b5


#set key samplen 4

#set key box vertical width 0.5 height 1 samplen 6 b l
set key box t r  
set key font ",35"
set key spacing 1.5

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "Length.eps"

p  [0:0.02]  "L6.txt" u (1/$1):($2) title "D=6"  with points  pointtype 2 lw 2  ps 8.0 lc rgb "#8B008B", "L7.txt" u (1/$1):($2) title "D=7"  with points  pointtype 4 lw 2  ps 8.0 lc rgb "#a40000","L9.txt" u (1/$1):($2) title "D=9"  with points  pointtype 8 lw 2  ps 8.0 lc rgb "#4e9a06","L12.txt" u (1/$1):($2+0.09) title "D=12"  with points  pointtype 10 lw 2  ps 8.0 lc rgb "#5c3566","L13.txt" u (1/$1):($2+0.09) title "D=13"  with points  pointtype 12 lw 2  ps 8.0 lc rgb "#ce5c00", x<0.015 ? G(x) : 1/0 notitle lw 4 lc rgb "#8B008B" dt 4, x<0.015 ? G1(x) : 1/0 notitle lw 4 lc rgb "#a40000" dt 4, x<0.015 ? G3(x) : 1/0 notitle lw 4 lc rgb "#8B008B" dt 4,x<0.015 ? G4(x) : 1/0 notitle lw 4 lc rgb "#5c3566" dt 4,x<0.012 ? G5(x) : 1/0 notitle lw 4 lc rgb "#ce5c00" dt 4  

#, x<0.015 ? G3(x) : 1/0 notitle lw 4 lc rgb "#4e9a06" dt 4, 

#"L8.txt" u (1/$1):($2) title "D=8"  with points  pointtype 6 lw 2  ps 8.0 lc rgb "#5c3566"

#, "L9.txt" u (1/$1):($2) title "D=9"  with points  pointtype 8 lw 2  ps 8.0 lc rgb "#4e9a06" ,"L12.txt" u (1/$1):($2) title "D=12"  with points  pointtype 10 lw 2  ps 8.0 lc rgb "#ff1493","L13.txt" u (1/$1):($2) title "D=13"  with points  pointtype 8 lw 2  ps 8.0 lc rgb "#c4a000"


# x<0.015 ? G2(x) : 1/0 notitle lw 4 lc rgb "#5c3566" dt 4


