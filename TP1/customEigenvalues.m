function ANS = customEigenvalues(A)

  tic;
  if (rows(A) == columns(A))
    H = customHassenbergDecomposition(A);
    ANS = customDoubleShiftQR(H);
  else
    printf ("Error: The given matrix is not squared.\n");
  endif
  toc;
endfunction
