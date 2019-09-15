set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "Schultz Iteration"   font ",40" textcolor rgb "black"
set ylabel "F_{sc}"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
set ytics nomirror ("10^{-1}" -1,"10^{-2}" -2,"10^{-3}" -3,"10^{-4}" -4,"10^{-5}" -5, "10^{-6}"-6, "10^{-7}" -7, "10^{-8}" -8, "10^{-9}" -9, "10^{-10}" -10, "10^{-11}" -11, "10^{-12}" -12,"10^{-13}" -13 );
set tics scale 4

set key box t l  
set key font ",40"
set key spacing 6


set output "FidelZAZ.eps"
p  "D3C3h2.txt" u ($1):(log(1-abs($2))/log(10))  t "D=3, {/Symbol c}=3" with points  pointtype 21 lw 3  ps 8.0 lc rgb "#c4a000","D3C9h2.txt" u ($1):(log(1-abs($2))/log(10))  t "D=3, {/Symbol c}=9" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#ef2929","D4C4h2.txt" u ($1):(log(1-abs($2))/log(10))  t "D=4, {/Symbol c}=4" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#ff1493","D4C10h2.txt" u ($1):(log(1-abs($2))/log(10))  t "D=4, {/Symbol c}=10" with points  pointtype 23 lw 3  ps 8.0 lc rgb "#ce5c00"


#, ,"D4AccRootX12.txt" u ($1):(log(1-abs($2)))  t "D=4, {/Symbol c}=12" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#c4a000"



#"D2AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=2, {/Symbol c}=4" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C"
#"D3AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=3, {/Symbol c}=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#8f5902"
#"D4AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=4, {/Symbol c}=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#5c3566"
