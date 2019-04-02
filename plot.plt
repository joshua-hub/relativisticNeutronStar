#!/usr/bin/gnuplot -p

plot "massvsradius.dat" using 2:3 title "Newtonian" w l, "massvsradius.dat" using 4:5 title "Relativistic" w l
set xlabel "Radius (km)"
set ylabel "Mass (solar Masses)"
set yrange [0:10]
set title "phys3071 as05 melsom 4259324"
set term postscript
set output "as05-42593249-melsom.ps"
replot
