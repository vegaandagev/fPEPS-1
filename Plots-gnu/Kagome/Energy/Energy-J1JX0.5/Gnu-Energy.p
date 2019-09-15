set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/D"   font ",35" textcolor rgb "black"
set ylabel "E"  font ",35"  textcolor rgb "black"

set xtics nomirror ("1/5" .2,"1/15" 0.066,"1/8" 0.125)
set tics scale 4

set key box t l  
set key font ",35"
set key spacing 6

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )


M(x)=a1+b1*((x**c1))
fit   M(x) "Echiral5.txt" u ($1):($2)    via a1,b1,c1



#M1(x)=a2+b2*(exp(-x*c2))
#fit   M1(x) "Echiral5.txt" u ($1):($2)    via a2,b2,c2

M2(x)=a3+b3*(x**(1))+c3*(x**(2))
fit   M2(x) "Echiral5.txt" u ($1):($2)    via a3,b3,c3

#G1(x)=0


set output "E.eps"
p [0.0:0.250]  [-0.48:-0.465] -0.4773 t "DMRG:L_y=12" with line lw 4 linetype 6 , "Echiral5.txt" u ($1):($2) t "simple-update" with points  pointtype 9 lw 3  ps 9.0 lc rgb "#DC143C", "Echiral5.txt"  u ($1):($6)  t "full-update"  with points  pointtype 8 lw 3  ps 9.0 lc rgb "#FF6347", x<0.335 ? M(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4,x<0.335 ? M2(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4


#x<0.27 ? M2(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4,   x<0.27 ? M(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4,x<0.27 ? M1(x) : 1/0 notitle lw 4 lc rgb "blue" dt 4
#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,
#"DMRG.txt" u ($1):($2) t "DMRG" with points pointtype 4 lw 3  ps 9.0 lc rgb "#4169E1"
