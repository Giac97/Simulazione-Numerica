set terminal pdf
set output 'comparison.pdf'

set multiplot layout 1,3 title "Initial and final comparison" font ",14"
set bmargin 5

set title "Best distance"
unset key
pl "distance.out" u 1 w l

set title "Initial"
unset key
pl "initial.out" u 1:2 w l

set title "Final"
unset key
plot "final.out" using 1:2 w l

unset multiplot

set output

