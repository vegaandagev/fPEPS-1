set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "{/Symbol l}"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
set ytics nomirror ("10^{-1}" -1,"10^{-2}" -2,"10^{-3}" -3,"10^{-4}" -4,"10^{-5}" -5, "10^{-6}"-6, "10^{-7}" -7, "10^{-8}" -8, "10^{-9}" -9, "10^{-10}" -10, "10^{-11}" -11, "10^{-12}" -12,"10^{-13}" -13 );
set tics scale 4

set key box b r  
set key font ",40"
set key spacing 6


set output "QRLinear.eps"
p [1:6] [-10:-0.6] "D2accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=2, l=2, 12{/Symbol \264}12" with points  pointtype 15 lw 3  ps 8.0 lc rgb "#ce5c00","l4D2accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=2, l=4, 12{/Symbol \264}12" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#a40000","l4D3accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=3, l=4, 12{/Symbol \264}12" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#ff1493","l4D2accIsingN8.txt" u ($1):(log(abs($2))/log(10))  t "D=2, l=4, 8{/Symbol \264}8" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#204a87", "D3accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=3, l=2, 12{/Symbol \264}12" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#ef2929"

#,"l4D2accIsingN8.txt" u ($1):(log(abs($2))/log(10))  t "D=2, l=4, 8{/Symbol \264}8" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#204a87","l4D3accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=3, l=4, 12{/Symbol \264}12" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#ce5c00"


#,"l4D2accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=2, l=4, 12{/Symbol \264}12" with points  pointtype 5 lw 3  ps 8.0 lc rgb "#f57900","D3accIsingN8.txt" u ($1):(log(abs($2))/log(10))  t "D=3, l=2, 8{/Symbol \264}8" with points  pointtype 2 lw 3  ps 8.0 lc rgb "#8f5902","l4D3accIsingN12.txt" u ($1):(log(abs($2))/log(10))  t "D=3, l=4, 12{/Symbol \264}12" with points  pointtype 12 lw 3  ps 8.0 lc rgb "#ce5c00"



#, "D3.txt" u ($1):(log(1-abs($3)))  t "D=3, l=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#C71585","D3.txt" u ($1):(log(1-abs($5)))  t "D=4, l=2" with points  pointtype 15 lw 3  ps 8.0 lc rgb "#204a87","D3.txt" u ($1):(log(1-abs($4)))  t "D=4, l=4" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#ce5c00"
