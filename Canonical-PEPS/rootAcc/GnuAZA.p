set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "Schultz Iteration"   font ",40" textcolor rgb "black"
set ylabel "F_{SC}"  font ",40"  textcolor rgb "black"

set tics scale 4

set key box t l  
set key font ",40"
set key spacing 1


set ytics nomirror ("10^{-1}" -1,"10^{-2}" -2,"10^{-3}" -3,"10^{-4}" -4,"10^{-5}" -5, "10^{-6}"-6, "10^{-7}" -7, "10^{-8}" -8, "10^{-9}" -9, "10^{-10}" -10, "10^{-11}" -11, "10^{-12}" -12,"10^{-13}" -13 );


set output "FidelZAZ.eps"
p   "D3AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=3, {/Symbol c}=4" with points  pointtype 15 lw 3  ps 4.0 lc rgb "#CD5C5C","D3AccRootX12.txt" u ($1):(log(1-abs($2)))  t "D=3, {/Symbol c}=12" with points  pointtype 17 lw 3  ps 4.0 lc rgb "#c4a000","D4AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=4, {/Symbol c}=4" with points  pointtype 19 lw 3  ps 4.0 lc rgb "#cc0000", "D4AccRootX12.txt" u ($1):(log(1-abs($2)))  t "D=4, {/Symbol c}=12" with points  pointtype 21 lw 3  ps 4.0 lc rgb "#ef2929",






