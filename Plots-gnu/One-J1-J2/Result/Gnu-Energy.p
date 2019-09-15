set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "J_{2}"   font ",40" textcolor rgb "black"
#set ylabel "E"  font ",40"  textcolor rgb "black"

#set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
#set xtics (0.6,0.61,0.62);
#set ytics (-0.471,-0.48,-0.475, -0.477,-0.478);


set key box t l  
set key font ",40"
set key spacing 8

#g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(1))
fit   G(x) "Energy.dat" u ($1):($3)    via a,b
#M(x)=a1+b1*(x)
#fit [0.10:0.167]  M(x) "J20.5.txt" u ($1):($2)    via a1,b1#,c
##M1(x)=a2+b2*(x**(1))+c2*(x**(2))
#fit   M1(x) "J20.5.txt" u ($1):($2)    via a2,b2,c2
#G1(x)=0


set output "HeisenbergJ2-5.eps"
p [0:0.35]  [-0.7:-0.38]  "Energy.dat"  title '2*2' with linespoints lt 4 dt 4 lw 3 pt 2 ps 8  lc rgb "#FF1493", "Energy.dat" u ($1):($3) title '3*2' with linespoints lt 4 dt 4 lw 3 pt 14 ps 8  lc rgb "#008080", "Energy.dat" u ($1):($4) title '4*3' with points lt 4 dt 4 lw 3 pt 10 ps 8  lc rgb "red", "Energy.dat" u ($1):($5) title '3*3' with linespoints lt 4 dt 4 lw 3 pt 10 ps 8  lc rgb "blue", G(x) notitle lw 4 lc rgb "#1E90FF" dt 4



