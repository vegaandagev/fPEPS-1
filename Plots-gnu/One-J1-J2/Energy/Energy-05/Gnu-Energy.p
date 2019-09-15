set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set x2tics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/L_{y}"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"
set x2label "1/D" font ",40" textcolor rgb "black"

set xtics nomirror ("1/6" .1666, "1/8" .125,"1/4" .25)
set x2tics ("1/6" .1666, "1/9" .111,"1/3" .333)

#set xtics (0.3,0.2,0.1);
#set ytics (-1.050,-1.484,-1.506,-1.55,-1.49,-1.515);
set tics scale 4
set ytics 0.02

set key box c r
set key font ",40"
set key spacing 6
set key width 2

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(2))+c*(x**(4))
fit [0.10:0.66]  G(x) "ipeps.txt" u ($1):($2)    via a,b,c
#M(x)=a1+b1*(x**(2))+c1*(x**(4))
#fit [0.10:0.66]  M(x) "ipeps.txt" u ($1):($4)    via a1,b1,c1
P(x)=a2+b2*(x**(3))#+c2*(x**(0))
fit [0.01:0.2]  P(x) "DMRG.txt" u ($1):($2)    via a2,b2#,c2



G1(x)=0


set output "HeisenbergJ2-5.eps"
p [0.0:0.450] [-1.5540:-1.480]   x<0.330 ? G(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4, "ipeps.txt" t "iPEPS" with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B", "DMRG.txt"  t "DMRG" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#FF1493", x<0.255 ? P(x) : 1/0 notitle lw 4 lc rgb "#C71585" dt 4


#"ipeps.txt"  u ($1):($4) t "iPEPS:J_2=0.548" with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8B008B"
#x<0.330 ? M(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4,
#"DMRG.txt"  u ($1):($3) t "DMRG" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#FF1493"
