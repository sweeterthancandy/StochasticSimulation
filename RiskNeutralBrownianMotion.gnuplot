set grid
set term pngcairo size 2560,1600

set datafile sep ','

set key autotitle columnhead

set macros 
POS = "at graph 0.92,0.9 font ',8'"

set output 'RiskNeutralBrownianMotion.png'
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set label 1 'a' @POS
plot \
        'RiskNeutralBrownianMotion.csv' u 1:3 w l , \
        for [i=8:27] 'RiskNeutralBrownianMotion.csv' u 1:i w l 
# --- GRAPH b
set label 1 'b' @POS
plot \
        for [i=8:22] 'RiskNeutralBrownianMotion.csv' u 1:(column(2)*column(i)) w l notitle
# --- GRAPH c
set label 1 'c' @POS
plot \
        'RiskNeutralBrownianMotion.csv' u 1:(10/column(2)) w l , \
        'RiskNeutralBrownianMotion.csv' u 1:3 w l , \
        'RiskNeutralBrownianMotion.csv' u 1:4 w l , \
        'RiskNeutralBrownianMotion.csv' u 1:5 w l , \
        'RiskNeutralBrownianMotion.csv' u 1:6 w l , \
        'RiskNeutralBrownianMotion.csv' u 1:7 w l , \
# --- GRAPH d
set label 1 'd' @POS
plot \
        'RiskNeutralBrownianMotion.csv' u 1:(10/column(2) - column(5)) w l title "Error_{1000}" , \
        'RiskNeutralBrownianMotion.csv' u 1:(10/column(2) - column(6)) w l title "Error_{2000}", \
        'RiskNeutralBrownianMotion.csv' u 1:(10/column(2) - column(7)) w l title "Error_{4000}"

unset multiplot

# vim:ft=gnuplot


