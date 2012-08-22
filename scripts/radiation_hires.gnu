set tics out
set xlabel "t"
set ylabel "$E$"
set xrange[0:10]
set yrange[0:0.8]
set term windows
set title ""
set arrow from 6.57,0.8 to 6.57,0 nohead lc rgb "black"
plot "../data/radiation_hires.dat" u ($1)/10:($3) linecolor rgb "black" pt 1 title "", "../data/radiation_hires_eff.dat" u ($1)/10:($3) linecolor rgb "black" pt 4 title ""
set term epslatex color
set output "radiation_hires.tex"
replot
set output