set grid
set term pngcairo size 2560,1600

set datafile sep ','

set key autotitle columnhead

set macros 
POS = "at graph 0.92,0.9 font ',8'"

set output 'BankAccount.png'
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set title "Interest rate realizations"
set label 1 'a' @POS
plot \
        for [i=9:11:2] 'BankAccount.csv' u 1:(column(i)) w l 
# --- GRAPH b
set title "Bank Account B(t)"
set label 1 'b' @POS
plot \
        for [i=10:20:2] 'BankAccount.csv' u 1:(column(i)) w l 
# --- GRAPH c
set label 1 'c' @POS
set title "Realizations of average value of call option, D(t)*V(t)"
plot \
        for [i=2:6] 'BankAccount.csv' u 1:(column(i)) w l 
# --- GRAPH d
set label 1 'd' @POS
set title "Discount Curve D(t) = 1/B(t)"
plot \
        for [i=4:6] 'BankAccount.csv' u 1:(1/column(i)) w l 

unset multiplot

# vim:ft=gnuplot


