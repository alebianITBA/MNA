function ANS = mvmran(n, nzr)

  t0 = time();
  ANS = zeros(n);
  position_in_column = [];
  rand_values = [];
  rand_values = rand(1, n * nzr);

  # Generate the row positions for each column
  for i = 1:n
    position_in_column = [ position_in_column, randperm(n, nzr) ];
  endfor

  # Now put each rand_values in the row and column that correspond
  index = 1;
  for i = 1:n
    for j = 1:nzr
      ANS(position_in_column(index), i) = rand_values(index);
      index = index + 1;
    endfor
  endfor
  printf("Time to create mvmran: %d\n", time() - t0);

endfunction
