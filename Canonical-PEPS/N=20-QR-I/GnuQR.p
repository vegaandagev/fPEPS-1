set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "Iteration"   font ",40" textcolor rgb "black"
set ylabel "F"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics nomirror (-0.4886,-0.4950, -0.4968,-0.4980);
set tics scale 4

set key box t r  
set key font ",40"
set key spacing 1


set ytics nomirror ("10^{-1}" -1,"10^{-2}" -2,"10^{-3}" -3,"10^{-4}" -4,"10^{-5}" -5, "10^{-6}"-6, "10^{-7}" -7, "10^{-8}" -8, "10^{-9}" -9, "10^{-10}" -10, "10^{-11}" -11, "10^{-12}" -12,"10^{-13}" -13 );


set output "QR.eps"
p [0:100]  "QRD3.txt" u ($1):(log(1-abs($2))/log(10))  t "D=3" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#DC143C", "QRD4.txt" u ($1):(log(1-abs($2))/log(10))  t "D=4" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#C71585", "QRD5.txt" u ($1):(log(1-abs($2))/log(10))  t "D=5" with points  pointtype 15 lw 3  ps 4.0 lc rgb "#fdc328"


#,"ACCD4QR.txt" u ($1):(log(1-abs($2)))  t "D=4, l=2" with points  pointtype 11 lw 3  ps 8.0 lc rgb "#ce5c00"




#"D2QR.txt" u ($1):(log(1-abs($2)))  t "D=2, l=2" with points  pointtype 6 lw 3  ps 8.0 lc rgb "#DC143C", "ACCD2QR.txt" u ($1):(log(1-abs($2)))  t "D=2, l=4" with points  pointtype 7 lw 3  ps 8.0 lc rgb "#C71585", 
