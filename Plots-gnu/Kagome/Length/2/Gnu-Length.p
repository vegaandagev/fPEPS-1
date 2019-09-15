set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/{/Symbol c}"   font ",35" textcolor rgb "black"
set ylabel "{/Symbol x}"  font ",40"  textcolor rgb "black"

#set xtics ("1/1000" 0.001,"1/250" 0.004, "1/100" 0.01)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);

#G(x)=b+a*x
#fit   G(x) "L6.txt" u (1/$1):($2)    via a,b

#G1(x)=b1+a1*x
#fit   G1(x) "L7.txt" u (1/$1):($2)    via a1, b1

G2(x)=b2+a2*x
fit   G2(x) "L8.txt" u (1/$1):($2)    via a2, b2

G3(x)=b3+a3*x
fit   G3(x) "L9.txt" u (1/$1):($2)    via a3, b3

G4(x)=b4+a4*x
fit  G4(x) "L10.txt" u (1/$1):($2)    via a4, b4

G5(x)=b5+a5*x
fit  G5(x) "L12.txt" u (1/$1):($2)    via a5, b5


G7(x)=b7+a7*x
fit  G7(x) "L13.txt" u (1/$1):($2)    via a7, b7


G8(x)=b8+a8*x
fit  G8(x) "L14.txt" u (1/$1):($2)    via a8, b8


G9(x)=b9+a9*x
fit  G9(x) "L15.txt" u (1/$1):($2)    via a9, b9


#set key samplen 4

#set key box vertical width 0.5 height 1 samplen 6 b l
set key box t r  
set key font ",35"
set key spacing 1.5

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "Length.eps"

p [0:0.02] [1.0:5.0]  "L8.txt" u (1/$1):($2) title "D=8"   with points  pointtype 6 lw 2  ps 8.0 lc rgb "#5c3566", "L9.txt" u (1/$1):($2) title "D=9"  with points  pointtype 8 lw 2  ps 8.0 lc rgb "#4e9a06", "L10.txt" u (1/$1):($2) title "D=10"  with points pointtype 14 lw 2  ps 8.0 lc rgb "#c4a000" ,"L12.txt" u (1/$1):($2) title "D=12"  with points  pointtype 10 lw 2  ps 8.0 lc rgb "#ff1493","L13.txt" u (1/$1):($2) title "D=13"  with points  pointtype 11 lw 2  ps 8.0 lc rgb "#FF1493","L14.txt" u (1/$1):($2) title "D=14"  with points  pointtype 14 lw 2  ps 8.0 lc rgb "#4682B4","L15.txt" u (1/$1):($2) title "D=15" with points  pointtype 16 lw 2  ps 8.0 lc rgb "#7FFFD4", x<0.015 ? G2(x) : 1/0 notitle lw 4 lc rgb "#5c3566" dt 4, x<0.015 ? G3(x) : 1/0 notitle lw 4 lc rgb "#4e9a06" dt 4, x<0.015 ? G4(x) : 1/0 notitle lw 4 lc rgb "#c4a000" dt 4,x<0.012 ? G7(x) : 1/0 notitle lw 4 lc rgb "#FF1493" dt 4,x<0.012 ? G8(x) : 1/0 notitle lw 4 lc rgb "#4682B4" dt 4,x<0.012 ? G9(x) : 1/0 notitle lw 4 lc rgb "#FF69B4" dt 4  ,x<0.012 ? G5(x) : 1/0 notitle lw 4 lc rgb "#00BFFF" dt 4 


