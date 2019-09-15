#set terminal postscript eps enhanced color font 'Helvetica,14' lw 3
set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,30' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}'  font ",40"  textcolor rgb "black"
set key b l 
set key box 6
set key font ",40"

#set key width 6

set key spacing 6


set output "Corrlog1.eps"

set xrange [log(2)/log(2):log(45)/log(2)]

set xtics ( 5, 15, 20,25,30,35, 40);
set ytics ( -3,-5,-7, -9,-11, -13,-15);


#M4(x)=(c4*x)+d4
#fit  [10:45] M4(x) "CorrelationH4.txt" u ($1):((log(abs($2)))*(1/log(2)))    via c4,d4 

M5(x)=(c5*x)+d5
fit  [20:45] M5(x) "CorrelationH5.txt" u ($1):((log(abs($2)))*(1/log(2)))    via c5,d5 

M6(x)=(c6*x)+d6
fit  [20:45] M6(x) "CorrelationH6.txt" u ($1):((log(abs($2)))*(1/log(2)))    via c6,d6 

M7(x)=(c7*x)+d7
fit  [20:45] M7(x) "CorrelationH7.txt" u ($1):((log(abs($2)))*(1/log(2)))    via c7,d7 

M8(x)=(c8*x)+d8
fit  [20:45] M8(x) "CorrelationH8.txt" u ($1):((log(abs($2)))*(1/log(2)))    via c8,d8 



Z5(x)=(a5*x)+b5
fit  [20:45] Z5(x) "CorrelationV5.txt" u ($1):((log(abs($2)))*(1/log(2)))    via a5,b5 

Z6(x)=(a6*x)+b6
fit  [20:45] Z6(x) "CorrelationV6.txt" u ($1):((log(abs($2)))*(1/log(2)))    via a6,b6 

Z7(x)=(a7*x)+b7
fit  [20:45] Z7(x) "CorrelationV7.txt" u ($1):((log(abs($2)))*(1/log(2)))    via a7,b7 

Z8(x)=(a8*x)+b8
fit  [20:45] Z8(x) "CorrelationV8.txt" u ($1):((log(abs($2)))*(1/log(2)))    via a8,b8 


#M8(x)=(c8*x)+d8
#fit  [10:40] M8(x) "CorrelationH8.txt" u ($1):((log(abs($2)))*(1/log(2)))    via c8,d8 




p  [5:40]  "CorrelationH5.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=5'  with points  pointtype 5 lw 3  ps 5.0 lc rgb "#FF69B4", "CorrelationH6.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationH7.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED","CorrelationH8.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=8' with points  pointtype 10 lw 3  ps 5.0 lc rgb "#6495ED","CorrelationV5.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=5'  with points  pointtype 5 lw 3  ps 5.0 lc rgb "red", "CorrelationV6.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=6'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "red","CorrelationV7.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=7' with points  pointtype 10 lw 3  ps 5.0 lc rgb "red","CorrelationV8.txt"  u ($1):((log(abs($2)))*(1/log(2))) t 'D=8' with points  pointtype 13 lw 3  ps 5.0 lc rgb "red", M5(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", M6(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",  M7(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", M8(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", Z5(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", Z6(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",  Z7(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",Z7(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "blue"

#,  M8(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970"
#, M4(x) notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970"


