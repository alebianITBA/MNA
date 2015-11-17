function cn = melcepstrum_transformation(frames, pr)
	for i = 1:columns(frames)
		cn(:,i) = melcepstrum(pr(:,i));
	endfor
endfunction