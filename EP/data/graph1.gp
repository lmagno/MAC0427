set xrange [0:2]
set yrange [0:2]
set xlabel "$x\_1$"
set ylabel "$x\_2$"

set term latex
set size square
set output "teste1.tex"

plot "teste1.dat" pt 7 ps 0.2, \
     1-x title "c1(x) = 0"

set output
