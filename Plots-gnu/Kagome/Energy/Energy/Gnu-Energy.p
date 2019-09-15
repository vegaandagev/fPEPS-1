set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "J_{chiral}"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key box b c  
set key font ",40"
set key spacing 6

#g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

#G(x)=a+b*(x**(2))+c*(x**(4))
#fit [0.10:0.66]  G(x) "J20.5.txt" u ($1):($2)    via a,b,c
#M(x)=a1+b1*(x)
#fit [0.10:0.167]  M(x) "J20.5.txt" u ($1):($2)    via a1,b1#,c
#M1(x)=a2+b2*(x**(1))+c2*(x**(2))
#fit   M1(x) "J20.5.txt" u ($1):($2)    via a2,b2,c2

G1(x)=0

set output "E.eps"

p [-0.2:0.2]    "Energy4.txt" u ($1):($2)  t "D=4" with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B", "Energy6.txt" u ($1):($2)  t "D=6" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#1E90FF","Energy8.txt" u ($1):($2)  t "D=8" with points  pointtype 4 lw 3  ps 8.0 lc rgb "#FF1493","Energy10.txt" u ($1):($2)  t "D=10" with points  pointtype 10 lw 3  ps 8.0 lc rgb "#2F4F4F"

#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,

