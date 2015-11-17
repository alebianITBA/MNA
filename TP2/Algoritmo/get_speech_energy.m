function energy = get_speech_energy(frame, frame_length)
energy =0;
	for n = 1:frame_length
		energy = energy + frame(n) ** 2;
	endfor
energy = log10(energy);
end