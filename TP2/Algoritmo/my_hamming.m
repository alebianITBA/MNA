function Ham = my_hamming(N)
	for i=1:N
		Ham(i) = 0.54 - 0.46*cos(2*pi*((i-1)/(N-1)));
	endfor
	Ham = Ham';
end