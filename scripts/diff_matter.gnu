set tics out
set xlabel "t"
set ylabel "$\\beta$"
set xrange[0:20]
set yrange[0:0.5]
set cbrange[0:0.05]
set palette defined (0 "white", 0.05 "#303030")
set term windows
set title "(b) matter"
plot "../data/diff_matter.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/(((1+x/32)**(0.66))**5+((1+x/32)**(0.66))**7) linecolor rgb "black" lw 6 title ""
set term epslatex color
set output "diff_matter.tex"
replot
set output
