set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/D"   font ",35" textcolor rgb "black"
set ylabel "m"  font ",35"  textcolor rgb "black"

set xtics nomirror ("1/5" .2,"1/16" 0.0625,"1/8" 0.125)
set tics scale 4

set key box t l  
set key font ",35"
set key spacing 6

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

#G(x)=a+b*(x**(1))+c*(x**(2))+d*(x**(3))
#fit G(x) "Echiral5.txt" u ($1):($4)    via a,b,c,d

M(x)=b1*((x**c1))
fit  M(x) "Echiral5.txt" u ($1):($4)    via b1,c1

G(x)=b*((x**c))
fit  G(x) "EKagome.txt" u ($1):($3)    via b,c

G2(x)=b2*((x**c2))
fit  G2(x) "Echiral5.txt" u ($1):($3)    via b2,c2


#[0.0:0.16]
#P(x)=exp(-x*d4)
#fit  P(x) "Echiral5.txt" u ($1):($3)    via d4


set output "M.eps"
p [0.0:0.250] [0:0.15]    "EKagome.txt"  u ($1):($7)  t "J_{ch}=0:full"  with points  pointtype 6 lw 3  ps 9.0 lc rgb "#9400D3", "EKagome.txt"  u ($1):($3)  t "J_{ch}=0:simple"  with points  pointtype 7 lw 3  ps 9.0 lc rgb "#cc0000","Echiral5.txt" u  ($1):($4) t "J_{ch}=0.5:full" with points  pointtype 8 lw 3  ps 9.0 lc rgb "#8b008b",  "Echiral5.txt" u  ($1):($3) t "J_{ch}=0.5:simple" with points  pointtype 9 lw 3  ps 9.0 lc rgb "#FF1493", x<0.27 ? M(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4, x<0.335 ? G(x) : 1/0 notitle lw 4 lc rgb "#cc0000" dt 4,x<0.335 ? G2(x) : 1/0 notitle lw 4 lc rgb "#cc0000" dt 4

#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,

