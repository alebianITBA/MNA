function frame_matrix_ret = frames(s,frame_length,frame_length_overlaped,Ham)

quantity_frames = floor((length(s) - frame_length) / (frame_length_overlaped + 1));
frame_matrix = ones(frame_length,quantity_frames-1);

frame_length_vec = [1:frame_length]';
step_vec = frame_length_overlaped .*[0:quantity_frames-1];
frame_matrix = step_vec + frame_length_vec;
frame_matrix_ret =s(frame_matrix).*Ham;

%To remove warnings, should use something like
%frame_matrix = bsxfun(@sum,step_vec , frame_length_vec);
%frame_matrix = bsxfun(@times, frame_matrix, Ham);
%we used warning ('off', 'Octave:broadcast'); in graphs() & fragmento2() & fragmento3() & fragmento4()

	%Equivalent code but more inefficient due to Octave's fors, we prefer to use vector sums
	%k=1;
	%for i =1:frame_length
	%	for j = 1:q
	%		frame_matrix(i,j) = k;
	%		k=k+frame_length_overlaped;
	%	endfor
	%	k=i+1;
	%endfor
	%size(frame_matrix) 160x234
end
