function pre_emphasized_speech = preemphasize(speech,a)
pre_emphasized_speech(1) = speech(1);
	for i=2:rows(speech)
		pre_emphasized_speech(i) = speech(i) - a * speech(i-1);
	endfor
end