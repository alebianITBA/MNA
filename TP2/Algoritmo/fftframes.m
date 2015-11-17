function fastfourierframes = fftframes(frames)
%Applying fft to each frame, each frame is represented by an entire column
for i=1:columns(frames)
	fastfourier(:,i)=fft(frames(:,i));
endfor
for i = 1:columns(frames)
	fastfourierframes(:,i) = (abs(fastfourier(:,i)).^2);
endfor
end