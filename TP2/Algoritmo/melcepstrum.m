function c = melcepstrum (s)
%Sk is the output power of the kth filter of the filter bank
%n is from 1 to 12

c = zeros(12,1);
n = rows(c);
z = zeros(rows(s),1);
K = rows(s);

for i= 1:n
	z(1)=log10(s(1))*cos(i*(1-0.5)*(pi/K));
	for k = 2:K
		z(k)= z(k-1) + log10(s(k))*cos(i*(k-0.5)*(pi/K));
	endfor	
	c(i) = z(K);
endfor
end