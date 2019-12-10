#!/usr/bin/gnuplot

##set term postscript enhanced monochrome lw 2 fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfss1200.pfb" "SFSS1200"
set term png enhanced
set out "plot_vs_reference.png"

set key bottom center
set xlabel "z / y"
set ylabel "U"

plot "sets/1/line_y_Ux.xy" u ($1-1):2 w l ls 2 t "OF", "sets/1/line_z_Ux.xy" u ($1-1):2 w l ls 2 notitle, "shercliff.dat" index 0 u 1:3 w l ls 1 title "Shercliff, 1953", "shercliff.dat" index 1 u 2:3 w l ls 1 notitle
