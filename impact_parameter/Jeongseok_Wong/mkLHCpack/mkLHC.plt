set term postscript eps enhanced color  # enhanced to get super & sub scripts
set xlabel "{/Symbol=20 D} {/Symbol=20 f}" 
set ylabel "{/Symbol=20 D} {/Symbol=20 h}" 
set zlabel " dN/d{/Symbol=20 D}{/Symbol=20 f}{/Symbol=20 D}{/Symbol=20 h}"   \
            offset 2 rotate by 90 
set xrange [ -2.5: 2.5 ]
set yrange [ -4.0: 4.0 ]
set zrange [ 0 : 0.02 ]
set ticslevel 0       
set hidden3d
set pm3d
set palette define (0. "blue", 0.003 "green",  0.008 "red")
#set palette rgbformulae 33,13,10
set title "pp at 7 GeV, 1 < p_T < 3 GeV/c" offset 0,-6 
splot "mkLHC.10" u 1:2:3 with lines notitle
set output "mkLHC10.ps"
replot
