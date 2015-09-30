set logscale y
set title "Residuals"
set ylabel 'Residual'
set xlabel 'Iteration'
#set term png
#set out "residuals.png"
plot "< cat log | grep 'Solving for Ux' | cut -d' ' -f9 | tr -d ','" title 'Ux' with lines,\
"< cat log | grep 'Solving for Uz' | cut -d' ' -f9 | tr -d ','" title 'Uz' with lines,\
"< cat log | grep 'Solving for pd' | cut -d' ' -f9 | tr -d ','" title 'p' with lines,\
"< cat log | grep 'Solving for epsilon' | cut -d' ' -f9 | tr -d ','" title 'epsilon' with lines,\
"< cat log | grep 'Solving for k' | cut -d' ' -f9 | tr -d ','" title 'k' with lines

