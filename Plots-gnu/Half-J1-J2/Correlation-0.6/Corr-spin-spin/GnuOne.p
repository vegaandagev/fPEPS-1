set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3
#test
set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "r"   font ",40" textcolor rgb "black"
set ylabel 'C^{s}, C^{d} '  font ",40"  textcolor rgb "black"




set key b l 
set key box 6
set key font ",40"
set tics scale 4
set format y ""
set format x ""
set logscale x 10
set logscale y 10
set ytics 1e-14, 10, 1
#set xtics 1e-10, 10, 1
set mytics 14
set mxtics 10
set ytics add ("10^{0}" 1.0,"10^{-2}" 0.01, "10^{-4}" 0.0001, "10^{-6}" 0.000001,"10^{-8}" 0.00000001, "10^{-10}" 0.0000000001, "10^{-12}" 0.000000000001)
set xtics add ("10^{0}" 1.0, "10^{-1}" 0.1, "10^{1}" 10,"10^{-2}" 0.01, "10^{-6}" 0.000001, "10^{-7}" 0.0000001,"10^{-8}" 0.00000001)

show logscale

set key spacing 6
set key spacing 10
set key width 2
set output "Corrlog.eps"

#M(x)=(c1*x)+d1
#fit  [log(25)/log(10):log(50)/log(10)] M(x) "CorrelationH7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    #via c1,d1 


#M2(x)=(c2*x)+d2
#fit  [log(30)/log(10):log(50)/log(10)] M2(x) "CorrelationV7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    #via c2,d2 



#M3(x)=(c3*x)+d3
#fit  [log(25)/log(10):log(45)/log(10)] M3(x) "CorrelationHH7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    #via c3,d3 


#M4(x)=(c4*x)+d4
#fit  [log(30)/log(10):log(45)/log(10)] M4(x) "CorrelationVV7.txt" u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10)))    #via c4,d4 



#G1(x)=0

p  [1.2:50] [0.0000000000001:1.0]   "CorrelationH8.txt"  u ($1):(abs($2)) t '{C^{s}}'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationV8.txt"  u ($1):(abs($2)) t '{C^{s}}'  with points  pointtype 4 lw 3  ps 5.0 lc rgb "#FF1493", "CorrelationHH.txt"  u ($1):(abs($2)) t '{C^{d}}' with points  pointtype 14 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationVV.txt"  u ($1):(abs($2)) t '{C^{d}}' with points  pointtype 14 lw 3  ps 5.0 lc rgb "#FF1493"

#,"CorrelationH7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) notitle  with points  pointtype 10 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationV7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) notitle  with points  pointtype 11 lw 3  ps 5.0 lc rgb "#FF1493", "CorrelationHH7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) notitle with points  pointtype 12 lw 3  ps 5.0 lc rgb "#3CB371","CorrelationVV7.txt"  u ((log(abs($1)))*(1/log(10))):((log(abs($2)))*(1/log(10))) notitle with points  pointtype 13 lw 3  ps 5.0 lc rgb "#FF1493"



#, x>1.4 ? M(x) : 1/0  notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970",x>1.4 ? M2(x) : #1/0 notitle with line dt 4 lw 4 linetype 4 lc rgb "#191970", x>1.3 ? M3(x) : 1/0  notitle #with line dt 4 lw 4 linetype 4 lc rgb "#191970",x>1.35 ? M4(x) : 1/0 notitle with line dt 4 #lw 4 linetype 4 lc rgb "#191970"



