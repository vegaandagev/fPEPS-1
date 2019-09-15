set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "iteration"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4

set key box t r 
set key font ",40"
set key spacing 6


set output "Energy.eps"

p [0:50] [-3.64:-3.51] "EnergyIter.txt" u ($1):($2)  t "D=2, {/Symbol c}=2" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#cc0000","EnergyIterImprove.txt" u ($1):($2)  t "D=2, {/Symbol c}=4" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#204a87", "EnergyIterD3.txt" u ($1):($2)  t "D=3, {/Symbol c}=3" with points  pointtype 15 lw 3  ps 8.0 lc rgb "#fdc328","EnergyIterD31.txt" u ($1):($2)  t "D=3, {/Symbol c}=6" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#a82296", -3.6140  t "DMRG" with line lw 4 linetype 7 lc rgb "#204a87"


#, -3.57  t "Matt" with line lw 4 linetype 7 lc rgb "#cd5c5c"
#"EnergyIterD3.txt" u ($1):($2)  t "D=3, {/Symbol g}=0.01" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#fdc328"


#"D3.txt" u ($1):(1-abs($2))  t "D=3, l=2" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C", "D3.txt" u ($1):(1-abs($3))  t "D=3, l=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#C71585",
