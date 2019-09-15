set terminal postscript eps size 6.3,6.3 enhanced color font 'Helvetica,25' lw 4

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "l_y"   font ",40" textcolor rgb "black"
#set ylabel "F"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
set ytics nomirror ( "2{/Symbol=30  \264}10^{-4}" 0.0002, "1{/Symbol=30  \264}10^{-4}" 0.0001, "1{/Symbol=30  \264}10^{-5}" 0.00001);
set tics scale 4
set boxwidth 0.6
set key box t l  
set key font ",40"
set key spacing 1


set output "QRLinear.eps"
#p  "D3.txt" u ($1):(0.5*log(1-abs($2)))  t "D=3, l=2" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#DC143C", "D3.txt" u ($1):(0.5*log(1-abs($3)))  t "D=3, l=4" with points  pointtype 17 lw 3  ps 4.0 lc rgb "#C71585","D3.txt" u ($1):(0.5*log(1-abs($5)))  t "D=4, l=2" with points  pointtype 15 lw 3  ps 4.0 lc rgb "#204a87","D3.txt" u ($1):(log(1-abs($4)))  t "D=4, l=4" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#f57900"

p [4:30]  "D2.txt" u ($1):((1-abs($2)))  t "D=2, l=2" with points  pointtype 11 lw 3  ps 5.0 lc rgb "#cc0000", "D2.txt" u ($1):((1-abs($3)))  t "D=2, l=4" with points  pointtype 9 lw 3  ps 5.0 lc rgb "#c71585", "D3.txt" u ($1):((1-abs($2))*0.5)  t "D=3, l=2" with points  pointtype 13 lw 3  ps 5.0 lc rgb "#f57900", "D3.txt" u ($1):((1-abs($3))) t "D=3, l=4" with points pointtype 15 lw 3  ps 5.0 lc rgb "#5cc863"


#"D3.txt" u ($1):(1-abs($2))  t "D=3, l=2" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C", "D3.txt" u ($1):(1-abs($3))  t "D=3, l=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#C71585",
