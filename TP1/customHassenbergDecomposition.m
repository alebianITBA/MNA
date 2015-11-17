function ANS = customHassenbergDecomposition(A)
  
  # Calculates the Hassenberg decomposition using Householder transforms
  ANS = A;
  n = rows(A);
  k = 1;
  for i = 1:(n - 2)
    x = ANS((k + 1):n, k);
    aux = n - (k+1) + 1;
    u = x - norm(x) * eye(aux)(:, 1);
    safe = [];
    if norm(u)^2 > 0
      safe = (2 * u * u') / norm(u)^2;
    else
      safe = 2 * u * u';
    endif
    p = eye(aux) - safe;
    aux = (k + 1) - 1;
    p = [ eye(n)(1:aux, :); [ zeros(n - 1)(aux:n-1, 1:aux), p ] ];
    k = k + 1;
    ANS = p * ANS * p;
  endfor
  # Round the elements of the matrix
  ANS = round(ANS .* 1000000) ./ 1000000;

endfunction
