set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "J_{c}"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol \266}^{2}H_{/Symbol c}/{/Symbol \266}J_{c}"  font ",40"  textcolor rgb "black"

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


set output "E-2D.eps"
p [.03:.32]   [-0.01:2] "Der7.txt" u ($1):(abs($3))  t "D=7" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C","Der9.txt" u ($1):(abs($3))  t "D=9" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#C71585","Der11.txt" u ($1):(abs($3))  t "D=11" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#8f5902"#,"DMRG.txt" u ($1):(abs($3))  t "DMRG" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#DA70D6"

#"Der10.txt" u ($1):(abs($3))  t "D=10" with points  pointtype 15 lw 3  ps 8.0 lc rgb "#8B008B"
#"Der8.txt" u ($1):(abs($3))  t "D=8" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#008B8B"
#"DMRGD.txt" u ($1):(abs($3))  t "DMRG" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#FF69B4"
#,"Der7.txt" u ($1):(abs($3))  t "D=7" with points  pointtype 10 lw 3  ps 8.0 lc rgb "#2F4F4F"
#"Der7.txt" u ($1):(abs($3))  t "D=7" with points  pointtype 5 lw 3  ps 8.0 lc rgb "#FF1493"
#-0.4958 t "Cluster PEPS" with line lw 4 linetype 6,

