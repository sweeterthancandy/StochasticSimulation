set grid
set term pngcairo size 2560,1600

set datafile sep ','

set key autotitle columnhead

set macros 
POS = "at graph 0.92,0.9 font ',8'"

set output 'CallOption.png'
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set title "Realizations of value of call option, V(t)"
set label 1 'a' @POS
plot \
        for [i=9:28] 'CallOption.csv' u 1:(column(i)) w l 
# --- GRAPH b
set title "Realizations of value of call option, D(t)*V(t)"
set label 1 'b' @POS
plot \
        for [i=9:28] 'CallOption.csv' u 1:(column(2)*column(i)) w l 
# --- GRAPH c
set label 1 'c' @POS
set title "Realizations of average value of call option, D(t)*V(t)"
plot \
        'CallOption.csv' u 1:3 w l, \
        for [i=4:8] 'CallOption.csv' u 1:(column(2)*column(i)) w l 
# --- GRAPH d
set label 1 'd' @POS
set title "Error(t) = BS(t) - D(T)*V(T)"
plot \
        'CallOption.csv' u 1:(column(3) - column(2)*column(6)) w l title "Error_{1000}" , \
        'CallOption.csv' u 1:(column(3) - column(2)*column(7)) w l title "Error_{2000}" , \
        'CallOption.csv' u 1:(column(3) - column(2)*column(8)) w l title "Error_{4000}" , \

unset multiplot

# vim:ft=gnuplot


