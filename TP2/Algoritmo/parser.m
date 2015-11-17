function ret = parser()
	fid = fopen('input.txt');
	if (fid==-1)
		ret={};
		return;
	endif


	str = fgetl(fid);
	if (!strcmp(str,'A:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if ((str(i)<'0' || str(i)>'9') && str(i) != '.')
			ret={};
			return;
		endif
	endfor
	ret{1} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'FS:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{2} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'BPS:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{3} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'FRAME LENGTH(ms):'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{4} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'overlappingcoef:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if ((str(i)<'0' || str(i)>'9') && str(i) != '.')
			ret={};
			return;
		endif
	endfor
	ret{5} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'lower_bank_filter_frec:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{6} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'higher_bank_filter_frec:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{7} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'quantity_of_banks:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{8} = sscanf(str,'%f',length(str));

	str = fgetl(fid);
	if (!strcmp(str,'nfft:'))
		ret={};
		return;
	endif
	str = fgetl(fid);
	for i=1:length(str)
		if (str(i)<'0' || str(i)>'9')
			ret={};
			return;
		endif
	endfor
	ret{9} = sscanf(str,'%f',length(str));


fclose(fid);

endfunction

