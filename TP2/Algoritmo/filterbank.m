function filter_banks = filterbank(lower_bank_filter_frec,higher_bank_filter_frec,quantity_of_banks,nfft,speech_rate)
%Mel transform of the highest frecuency(Hz)-> Mel
high = mel_mapping(higher_bank_filter_frec/2);
%Mel transform of the lowest frecuency(Hz)-> Mel
low = mel_mapping(lower_bank_filter_frec);
%Range between highest and lowest in Mel scale
dif = (high-low)/(quantity_of_banks-1);

%N points -> N-2 filterbanks, in Mel frecuency
for i= 1:quantity_of_banks
	vec(i) =low + (i-1)*dif;
endfor

%N points -> N-2 filterbanks, in Hz frecuency
for j=1:quantity_of_banks
	anti_vec(j) = mel_inverse(vec(j));
endfor

%Round those frequencies to the nearest FFT bin
for k=1:quantity_of_banks
	f(k) = floor((nfft+1)*anti_vec(k)/speech_rate);
endfor


%(Quantity of banks - 2) as 2 points are the extreme ones, the high and low
filter_banks = zeros(quantity_of_banks-2,nfft/2+1);

for m = 2:quantity_of_banks-2	
from = f(m-1);
to = f(m); 
	for k = from:to
		filter_banks(m,k) = (k - f(m-1))/(f(m)-f(m-1));
	end
endfor

for m = 1:quantity_of_banks-2
from = f(m);
to = f(m+1);
	for k = from:to
		filter_banks(m,k) = (f(m+1)-k)/(f(m+1)-f(m));
	end
endfor

end