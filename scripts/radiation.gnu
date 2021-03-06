set tics out
set xlabel "t"
set ylabel "$\\beta$"
set xrange[0:20]
set yrange[0:0.5]
set cbrange[0:1]
set palette defined (0 "white", 1 "#303030")
set term windows
set title "(a) radiation"
plot "../data/radiation.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/(sqrt(1+x/32)**5+sqrt(1+x/32)**7) linecolor rgb "black" lw 6 title ""
set term epslatex color
set output "radiation.tex"
replot
set output
