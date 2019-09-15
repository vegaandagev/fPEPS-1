set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "position x"   font ",40" textcolor rgb "black"
#set ylabel "{/Symbol D}"  font ",40"  textcolor rgb "black"

#set xtics nomirror ("1/5" .2, "1/9" 0.11111,"1/3" .333)
#set xtics (0.3,0.2,0.1);
set ytics nomirror ("10^{-1}" -1,"10^{-2}" -2,"10^{-3}" -3,"10^{-4}" -4,"10^{-5}" -5, "10^{-6}"-6, "10^{-7}" -7, "10^{-8}" -8, "10^{-9}" -9, "10^{-10}" -10, "10^{-11}" -11, "10^{-12}" -12,"10^{-13}" -13 );
set tics scale 4

set key box b c  
set key font ",40"
set key spacing 1

set key at 0, -6
set output "QR.eps"
p [-10:10]   "D3NormN16.txt" u ($1-8):(log(abs($2-$3)*0.006/($2))/log(10))  t "{/Symbol l}=2.00, l=2" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#a40000","D3NormN16.txt" u ($1-8):(log(abs($2-$3)*0.5/($2))/log(10))  t "{/Symbol l}=3.05, l=2" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#5c3566","h3D3NormN16.txt" u ($1-8):(log(abs($2-$3)/($2))/log(10))  t "{/Symbol l}=3.05, l=4" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ce5c00","h3D3NormN16.txt" u ($1-8):(log(abs($2-$3)*0.0007/($2))/log(10))  t "{/Symbol l}=2.00, l=4" with points  pointtype 15 lw 3  ps 4.0 lc rgb "#cd5c5c"



#,"ACCD2NormN12.txt" u ($1-6):(log(abs($2-$3)/($2)))  t "D=2, l=4" with points  pointtype 10 lw 3  ps 8.0 lc rgb "#cc0000","D2NormN12.txt" u ($1-6):(log(abs($2-$3)/($2)))  t "D=2, l=2" with points  pointtype 12 lw 3  ps 8.0 lc rgb "#cc0000"





#,"ACCD3NormN12.txt" u ($1-6):(log(abs($2-$3)/($2)))  t "D=3, l=4" with points  pointtype 15 lw 3  ps 8.0 lc rgb "#cc0000"


#"D2NormN16.txt" u ($1-8):(log(abs($2-$3)/($2)))  t "D=2, l=2" with points  pointtype 13 lw 3  ps 8.0 lc rgb "#C71585",


#,"D3.txt" u ($1):(log(1-abs($2)))  t "D=3, l=2" with points  pointtype 9 lw 3  ps 8.0 lc rgb "#DC143C", "D3.txt" u ($1):(log(1-abs($3)))  t "D=3, l=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#C71585","D3.txt" u ($1):(log(1-abs($5)))  t "D=4, l=2" with points  pointtype 15 lw 3  ps 8.0 lc rgb "#204a87", "D3.txt" u ($1):(log(1-abs($6)))  t "D=4, l=4" with points  pointtype 17 lw 3  ps 8.0 lc rgb "#c17d11"  
