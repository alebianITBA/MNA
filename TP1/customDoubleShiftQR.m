function ANS = customDoubleShiftQR(H)

  max_iterations = 50;
  convergence = 0.0001;
  eigenvalues = [];
  rows = rows(H);

  do
    for i = 1:max_iterations
      [Q,R] = customQR(H);
      H = R*Q;
      if abs(H(rows, rows-1)) < convergence
        eigenvalues(end+1) = H(rows, rows);
        H = H(1:rows-1, 1:rows-1);
        rows = rows - 1;
        break;
      endif

      if i == max_iterations
        submatrix = H(rows-1:rows, rows-1:rows);
        b = -1 * (submatrix(1, 1) + submatrix(2, 2));
        c = (submatrix(1, 1) * submatrix(2, 2)) - (submatrix(1, 2) * submatrix(2, 1));

        determ = b^2 - 4*c;
        if determ >= 0
          eigenvalues(end+1) = (-1*b + sqrt(determ)) / 2;
          eigenvalues(end+1) = (-1*b - sqrt(determ)) / 2;
        end

        if determ < 0
          eigenvalues(end+1) = complex(-1*b/2, sqrt(determ*-1) / 2);
          eigenvalues(end+1) = complex(-1*b/2, -1*sqrt(determ*-1) / 2);
        end

        H = H(1:rows-2, 1:rows-2);
        rows = rows - 2;
      endif
    endfor
  until rows < 2

  eigenvalues(end+1) = H;
  ANS = eigenvalues;

endfunction
