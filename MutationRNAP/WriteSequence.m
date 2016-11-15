function WriteSequence(path, Sequence)
	fout = fopen(path, 'w');
	fprintf(fout, '%s\n', Sequence);
	fclose(fout);
end