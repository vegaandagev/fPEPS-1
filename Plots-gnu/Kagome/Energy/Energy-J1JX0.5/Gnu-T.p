set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/D"   font ",35" textcolor rgb "black"
set ylabel "{/Symbol D}T"  font ",35"  textcolor rgb "black"

set xtics nomirror ("1/5" .2,"1/16" 0.0625,"1/8" 0.125)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4
#-0.49580,

set key box t l  
set key font ",35"
set key spacing 7
set key width 2

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

P(x)=b2*(x**c2)
fit P(x) "Delta.txt" u ($1):($7)    via b2,c2

M(x)=b1*((x**c1))
fit  M(x) "Delta.txt" u ($1):($5)    via b1,c1

#G(x)=b*((x**c))
#fit  G(x) "Delta.txt" u ($1):($3)    via b,c


#[0.0:0.16]
#P(x)=exp(-x*d4)
#fit  P(x) "Echiral5.txt" u ($1):($3)    via d4


set output "Delta.eps"
p  [0.0:0.250]  "Delta.txt"  u ($1):($5)  t  "{/Symbol D}T_{x}:J_{ch}=0.5"  with points  pointtype 9 lw 3  ps 9.0 lc rgb "#ff1493", "Delta.txt" u  ($1):($7) t "{/Symbol D}T_{x-y}:J_{ch}=0.5" with points pointtype 8 lw 3  ps 9.0 lc rgb "#3465a4", x<0.27 ? M(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4, x<0.27 ? P(x) : 1/0 notitle lw 4 lc rgb "#cc0000" dt 4

#, x<0.27 ? P(x) : 1/0 notitle lw 4 lc rgb "#729fcf" dt 4


#"Delta.txt" u  ($1):($6) t "{/Symbol D}T_{y}:J_{{/Symbol c}}=0.5" with points  pointtype 4 lw 3  ps 8.0 lc rgb "#8b008b"



