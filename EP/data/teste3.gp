xmin = -1
xmax = 3
ymin = -1
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

c1(x, y) = x - 1
c2(x, y) = y

set table '../aux/c1_3'
splot c1(x, y)
unset table

set table '../aux/c2_3'
splot c2(x, y)
unset table

set term epslatex
set output "../aux/teste3.tex"

set xzeroaxis lt -1 lw 3
set xtics axis (1) offset 0.5,0
set label "$x$" at xmax-0.1,-.1

set yzeroaxis lt -1 lw 3
unset ytics
set label "$y$" at -.1,ymax-0.1 right


set size square
plot "../aux/teste3.dat" pt 1 ps 1 t "", \
     "../aux/c1_3" u 1:2 title "$x = 1$" w lines lw 3,\
     "../aux/c2_3" u 1:2 title "$y = 0$" w lines lw 3
