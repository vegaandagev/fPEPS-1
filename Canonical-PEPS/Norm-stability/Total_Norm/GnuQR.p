set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "{/Symbol a}"   font ",40" textcolor rgb "black"
set ylabel "N"  font ",40"  textcolor rgb "black"
set xtics (0.01,0.02,0.03,0.04,0.05);

set tics scale 4

set key box t c  
set key font ",35"
set key spacing 1

set output "Random.eps"
#p [0.0:0.03]  "NormD2R.txt" u ($1):($3)  t "random PEPS" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#ef2929","NormD2R.txt" u ($1):($2)  t "random cPEPS" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ce5c00"


#, "NormD3L35.txt" u ($1):($2)  t "random PEPS" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#f57900"


p "NormD3L25.txt" u ($1):($3)  t "PEPS, {/Symbol l}=0.25" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#ef2929","NormD3L25.txt" u ($1):($2)  t "cPEPS, {/Symbol l}=0.25" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ce5c00", "NormD3L35.txt" u ($1):($3)  t "PEPS, {/Symbol l}=3.5" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#5c3566", "NormD3L35.txt" u ($1):($2)  t "cPEPS, {/Symbol l}=3.5" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#f57900"

#,"NormD4L25.txt" u ($1):($2)  t "cPEPS, D4" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#ef2929", "NormD4L25.txt" u ($1):($3)  t "PEPS, D4" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#5c3566", "NormD3L25.txt" u ($1):($2)  t "cPEPSA, D3" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#f57900"
