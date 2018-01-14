set term postscript color enhanced eps size 3,2.4 font 'Helvetica-Italic, 14'
set output '../figures/run003_initial_data.eps'
set xlabel 'X(u)'
set ylabel 'Y(u)'

plot[-2:2][-2:2] \
	'../../data/test.disc.file' u 2:3 t 'two poles initial data' w l

set output

