set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "{/Symbol c}"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"

#set xtics nomirror (20,40,60,80)

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics nomirror (100,200,300,400)

set tics scale 4

set key box t r 
set key font ",35"
#set key spacing 6


set output "E.eps"


p  "ED4.txt" u  ($1):(log((abs((($2)+10.642663)/10.642663)))/log(10))  t "D=4, V=0.1"  with  points  pointtype  5  lw 3  ps 4.0 lc rgb "#cc0000","ED4.txt" u  ($1):(log((abs((($3)+10.72835)/10.72835)))/log(10))  t "D=6, V=0.1"  with  points  pointtype  7  lw 3  ps 4.0 lc rgb "#ce5c00","ED4.txt" u  ($1):(log((abs((($4)+0.50535114070)/0.50535114070)))/log(10))  t "D=6, V=2.0"  with  points  pointtype  9  lw 3  ps 4.0 lc rgb "#a40000"


#,"ED4.txt" u  ($1):(log((abs((($3)+0.6651625)/0.6651625)))/log(10))  t "D=4, V=0.1"  with  points  pointtype  7  lw 3  ps 4.0 lc rgb "#ce5c00","ED5.txt" u  ($1):(log((abs((($2)+0.663703025931)/0.663703025931)))/log(10))  t "D=5, V=0.1, U(1), D"  with  points  pointtype  9  lw 3  ps 4.0 lc rgb "#a40000","ED5.txt" u  ($1):(log((abs((($3)+0.663703025931)/0.663703025931)))/log(10))  t "D=5, V=0.1, U(1), S"  with  points  pointtype  11  lw 3  ps 4.0 lc rgb "#c4a000"





