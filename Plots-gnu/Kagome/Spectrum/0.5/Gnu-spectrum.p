set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,20' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "{/Symbol f}"   font ",40" textcolor rgb "black"
#set ylabel "-{/Symbol l}"  font ",40"  textcolor rgb "black"

#set xtics ("1/5" .2, "1/8" .125,"1/3" .333)
#set xtics (0.3,0.2,0.1);
#set ytics (-0.4886,-0.4950,-0.49580, -0.4968,-0.4980);


set key font ",35"
set key spacing 6

set key box b l  

#set style fill transparent solid 0.5 noborder
#set style circle radius 0.02
set output "spectrum.eps"
p  [-2.150:2.150] [0:1.5]  "SpectrumD15.txt"  u ($3 + 0.6):4  title  "D15"   with points  pointtype 9 lw 2 ps 8.0 lc rgb "#FF0000","SpectrumD12.txt"  u ($3 + 0.2):4  title  "D12"  with  points  pointtype 4 lw 2 ps 8.0 lc rgb "#6A5ACD"

#, ,"SpectrumD10.txt"  u ($3 + 0.07):4  title  "D10"   with points  pointtype 4 lw 2 ps 8.0 lc rgb "#6A5ACD"


#"SpectrumD13.txt"  u 3:4  title  "D13"  with points  pointtype 2 lw 2 ps 6.0 lc rgb "#4169E1","SpectrumD6.txt"  u 3:4  title  "D6"   with points  pointtype 4 lw 2 ps 6.0 lc rgb "#FF1493","SpectrumD8.txt"  u 3:4  title  "D8"   with points  pointtype 6 lw 2 ps 6.0 lc rgb "#FF69B4",
