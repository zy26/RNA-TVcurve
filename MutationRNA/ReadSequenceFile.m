function seq = ReadSequenceFile(input)
	fid = fopen(input);
	
	tline = fgetl(fid);
	while ischar(tline)
		disp(tline);
		seq = tline;
		tline = fgetl(fid);
	end

	fclose(fid);

	seq = upper(seq);
	seq = strrep(seq, 'T', 'U');
end