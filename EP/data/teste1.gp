set term epslatex
set output "teste1.tex"

min = -1
max = 2
set xrange [min:max]
set yrange [min:max]

unset border

set xzeroaxis lt -1 lw 3
set xtics axis (1)
set label "$x$" at max-0.1,-.1

set yzeroaxis lt -1 lw 3
set ytics axis (1)
set label "$y$" at -.1,max-0.1


set label "$(\\frac{1}{2},\\frac{1}{2})$" at 0.55,0.55 left
set size square
plot "teste1.dat" pt 7 ps 0.5 title "", \
     1-x w lines lw 3 title "$x + y = 1$"
