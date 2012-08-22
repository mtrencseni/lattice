set tics out
set xrange[0:20]
set yrange[0:0.5]
set xlabel "t"
set ylabel "$\\beta$"
set cbrange[0:0.05]
set palette defined (0 "white", 0.05 "#303030")
set term windows
set title "(c) lambda"
plot "../data/diff_lambda.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/((exp(1+x/32)/exp(1))**5+(exp(1+x/32)/exp(1))**7) lc rgb "black" lw 6 title ""
set term epslatex color
set output "diff_lambda.tex"
replot
set output
