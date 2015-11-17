function deltas_answer = deltas(cn)
	deltas_answer(1) = (2*(cn(3)-cn(1)) + (cn(2)-cn(1))) / 10;
	deltas_answer(2) = (2*(cn(4)-cn(1)) + (cn(3)-cn(1))) / 10;
	for t = 3:11
		deltas_answer(t) = (2*(cn(t+2)-cn(t-2)) + (cn(t+1)-cn(t-1))) / 10;
	endfor
	deltas_answer(12) = (2*(cn(13)-cn(10)) + (cn(13)-cn(11))) / 10;
	deltas_answer(13) = (2*(cn(13)-cn(11)) + (cn(13)-cn(12))) / 10;
	
endfunction