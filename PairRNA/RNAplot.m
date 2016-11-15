function RNAplot(sequence, stru, output)

RNAplotFile = tempname;
fout = fopen(RNAplotFile, 'w');
fprintf(fout, '%s\n', sequence);
fprintf(fout, '%s\n', stru);
fprintf(fout, '%s\n', '@');
fclose(fout);

RunRNAplotCmd('RNAplot -o svg', RNAplotFile)

movefile('rna.svg', strcat(output, '.svg'));

delete(RNAplotFile);