set cntrparam levels discrete 0
set view map
unset surface
set isosamples 1000,1000
set xrange [-3:4]
set yrange [-2:2]
set xlabel "x"
set ylabel "y"
set contour

c1(x, y) = (x-1)**2 + (y-1)**2 - 1
c2(x, y) = (x-1)**2 + (y+1)**2 - 1

set table 'c1_2'
splot c1(x, y)
unset table 

set table 'c2_2'
splot c2(x, y)
unset table 


set term pdfcairo font "inconsolata,12"
set size ratio 0.5
set output "teste2.pdf"

set arrow from -3,0 to 4,0 nohead
set arrow from 0,-2 to 0,2 nohead
plot "teste2.dat" pt 7 ps 0.2, \
     "c1_2" u 1:2 title "c1(x) = 0" w lines,\
     "c2_2" u 1:2 title "c2(x) = 0" w lines
     
set output