set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/D"   font ",40" textcolor rgb "black"
#set ylabel "m"  font ",40"  textcolor rgb "black"

set xtics ("1/5" .2, "1/9" .1111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.473,-0.477,-0.4782, -0.482);


set key box b r  
set key font ",40"
set key spacing 6

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(1))
fit [0.10:0.31]  G(x) "ipeps5.txt" u ($1):($3)    via a,b
M(x)=a1+b1*(x**(1))
fit [0.10:0.2]  M(x) "ipeps55.txt" u ($1):($5)    via a1,b1
C(x)=a2+b2*(x**(1))
fit [0.10:0.31]  C(x) "ipeps545.txt" u ($1):($5)    via a2,b2


set output "Mag.eps"
p [0.0:0.40]    "ipeps5.txt" u ($1):($3)  t 'J_2=0.5' with points  pointtype 2 lw 3  ps 8.0 lc rgb "#DC143C","ipeps545.txt" u ($1):($5)  t 'J_2=0.545' with points  pointtype 4 lw 3  ps 8.0 lc rgb "#FF6347","ipeps55.txt" u ($1):($5)  t 'J_2=0.55' with points  pointtype 8 lw 3  ps 8.0 lc rgb "#DC143C","ipeps548.txt" u ($1):($6)  t 'J_2=0.548' with points  pointtype 14 lw 3  ps 8.0 lc rgb "#FF6347", x<0.330 ? M(x) : 1/0 notitle lw 4 lc rgb "#FF00FF" dt 4, x<0.330 ? C(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4, x<0.330 ? G(x) : 1/0 notitle lw 4 lc rgb "#FF00FF" dt 4


#"ipeps55.txt" u ($1):($4)  t 'J_2=0.55' with points  pointtype 6 lw 3  ps 8.0 lc rgb "#BA55D3",
