xmin = -1
xmax = 3
ymin = -2
ymax = 2
set xrange [xmin:xmax]
set yrange [ymin:ymax]

set key tmargin
unset border

set cntrparam levels discrete 0
set view map
unset surface
set isosamples 1000,1000
set contour

c1(x, y) = (x-1)**2 + (y-1)**2 - 1
c2(x, y) = (x-1)**2 + (y+1)**2 - 1

set table '../aux/c1_2'
splot c1(x, y)
unset table

set table '../aux/c2_2'
splot c2(x, y)
unset table

set term epslatex
set output "teste2.tex"

set xzeroaxis lt -1 lw 3
set xtics axis (1)
set label "$x$" at xmax-0.1,-.1

set yzeroaxis lt -1 lw 3
unset ytics
set label "$y$" at -.1,ymax-0.1 right


set size square
plot "teste2.dat" pt 1 ps 1 t "", \
     "../aux/c1_2" u 1:2 title "$(x-1)^2+(y-1)^2 = 1$" w lines lw 3,\
     "../aux/c2_2" u 1:2 title "$(x-1)^2+(y+1)^2 = 1$" w lines lw 3
