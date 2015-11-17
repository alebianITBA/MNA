function [Q, R] = customQR(H)

  n = rows(H);
  R = H;
  Q = eye(n);
  for i = 1:n
    for j = i+1:n
      if R(j,i) == 0
        c = 1;
        s = 0;
      elseif abs(R(j,i)) < abs(R(i,i))
        t = R(j,i) / R(i,i);
        c = 1 / sqrt(1 + t^2);
        s = c*t;
      else
        z = R(i,i) / R(j,i);
        s = 1 / sqrt(1 + z^2);
        c = s*z;
      endif
      G = [c,s;-1*s,c];
      R([i,j],i:n) = G * R([i,j],i:n);
      Q([i,j],1:n) = G * Q([i,j],1:n);
    endfor
  endfor
  Q = Q';

endfunction
