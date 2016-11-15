function SaveStructure(header, sequence, stru, output)

RNAplotFile = strcat(output, '.txt');
fout = fopen(RNAplotFile, 'w');
fprintf(fout, '> %s\n', header);
fprintf(fout, '%s\n', sequence);
fprintf(fout, '%s\n', stru);
fclose(fout);