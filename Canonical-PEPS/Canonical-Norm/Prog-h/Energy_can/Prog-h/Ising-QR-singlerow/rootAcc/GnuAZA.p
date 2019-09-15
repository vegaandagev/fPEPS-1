set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "Schultz Iteration"   font ",40" textcolor rgb "black"
set ylabel "F_{Schultz}"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4

set key box t l  
set key font ",40"
set key spacing 6


set output "FidelZAZ.eps"
p  "D2AccRootX12.txt" u ($1):(log(1-abs($2)))  t "D=2, {/Symbol c}=12" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#C71585", "D3AccRootX12.txt" u ($1):(log(1-abs($2)))  t "D=3, {/Symbol c}=12" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#CD5C5C","D4AccRootX12.txt" u ($1):(log(1-abs($2)))  t "D=4, {/Symbol c}=12" with points  pointtype 19 lw 3  ps 8.0 lc rgb "#c4a000"



#"D2AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=2, {/Symbol c}=4" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C"
#"D3AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=3, {/Symbol c}=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#8f5902"
#"D4AccRootX4.txt" u ($1):(log(1-abs($2)))  t "D=4, {/Symbol c}=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#5c3566"
