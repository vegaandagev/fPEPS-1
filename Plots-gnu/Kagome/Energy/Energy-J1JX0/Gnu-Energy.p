set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/D"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"

set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key b r 
set key box 4
set key font ",38"
set key spacing 8
set output "Length.eps"
set tics scale 3
set tics scale 3
set key width 4

#P(x)=(x**(v1))*z1+b1
M(x)=b1*(x**(a1))+c1
#M(x)=exp(x*a1)#+b1
#M(x)=exp(x*a1)


#fit  M(x) "EKagome.txt"  u ($1):($2)  via c1,d1 ,b1
fit  M(x) "EKagome.txt"  u ($1):($2)  via a1, b1,c1 
#fit  P(x) "EKagome.txt"  u ($1):($2)  via v1, z1#,c1 
#fit  M(x) "EKagome.txt"  u ($1):($2)  via a1#, b1#,c1 

G(x)=(x**(c2))*d2+m1
#G(x)=exp(x*d2)#+c2
#G(x)=exp(x*d2)+c2
fit  G(x) "Echiral5.txt"  u ($1):($2)  via c2,d2,m1
#fit  G(x) "Echiral5.txt"  u ($1):($2)  via d2


set output "Energy.eps"
p [0.0:0.50] [-0.5:-0.4]   "EKagome.txt"  u ($1):($2)  t "J_{{/Symbol c}}=0"  with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B","Echiral5.txt"  u ($1):($2)  t "J_{{/Symbol c}}=0.5"  with points  pointtype 4 lw 3  ps 8.0 lc rgb "#ce5c00", M(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2", G(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#8A2BE2", G(x) notitle with line dt 3 lw 3 linetype 3 lc rgb "#ce5c00"

#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,

