%[speech,fs,bps] = wavread("daniel_fragmento_1.wav");
function melcepstrum = mfcc(speech, bps)
	
	input_user = parser();

	if (length(input_user)==0)
		printf('Error en los parametros. Revisar input.txt\n');
		return;
	endif


	%VARIABLES SECTION
		%pre-emphasize coefficient
		a = input_user{1}; 
		%sample rate in hz
		fs = input_user{2} ;
		%number of bits per sample
		bps = input_user{3}; 
		%frame length in ms
		frame_length_ms = input_user{4}; 
		%coefficient of overlap between 2 consecutive frames
		overlapping_coef = input_user{5};
		
		%Lower frecuency used in Mel filterBank
		lower_bank_filter_frec = input_user{6};
		%Higher frecuency used in Mel filterBank
		higher_bank_filter_frec = input_user{7};
		%160 samples as it was recorded at 8000hz with a 
		%frame_length_ms = 0,02s = 20ms to guarantee stationary
		frame_length = fs * (frame_length_ms/1000);
		%80 samples as it was recorded at 8000hz with a 
		%frame_length_ms = 0,01s = 10ms of overlap
		frame_length_overlaped = fs * overlapping_coef * (frame_length_ms/1000);
		%quantity of filter banks
		quantity_of_banks = input_user{8};

		nfft = input_user{9};

	%ALGORITHM SECTION:

		%PRE-EMPHASIS BLOCK: emphasizing the signal
			pre_emphasized_speech = preemphasize(speech,a);
		
		%WINDOW AND OVERLAP BLOCK:
			%Rename variable to be consistent with paper											
			N = frame_length;
			%Own implementation of Hamming-function to reduce high 
			%frecuency noise introduced by framing process
			Ham = my_hamming(N);
			%Once emphasized, it is divided into frames
			frames = frames(pre_emphasized_speech,frame_length,frame_length_overlaped,Ham);
	
		% %FFT BLOCK: Apply fft to each frame 
		 	fastfourier = fftframes(frames);
	
		% %MEL FRECUENCY AND FILTER BANK BLOCK:
		 	filtered = filterbank(lower_bank_filter_frec,higher_bank_filter_frec,quantity_of_banks,rows(frames),fs);
		
		% %CEPSTRUM BLOCK:
			C = melcepstrum_transformation(frames, filtered * fastfourier(1:frame_length_overlaped+1, :));

		% %LOGGED ENERGY BLOCK: 	
		 	for i=1:columns(frames)
		 		C(13,i) = get_speech_energy(frames(:,i),frame_length);
		 	endfor
		
		% %DELTA BLOCK:
			for i=1:columns(frames)
			 	deltass = deltas(C(:,i));
			 	for j = 1:13
					C(j+13,i) = deltass(j);
				endfor
			endfor
			
			
		melcepstrum = C;

end