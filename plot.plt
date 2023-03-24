unset logscale x
set xlabel "X" 
set ylabel "Y"
set grid
set yrange [-0.5:1]
set pm3d
set terminal png enhanced truecolor 
  set output 'tempfile.png'            
    plot  "./dat2.txt" with lines title "МНРП", "./dat1.txt" with lines title "Исходная функция"
  set out                              # restore the output redirection
set terminal GNUTERM                    