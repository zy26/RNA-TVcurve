function [RNApdist] = RNApdist(sequence1, sequence2)

RNApdistPair = tempname;
fout = fopen(RNApdistPair, 'w');
fprintf(fout, '%s\n', sequence1);
fprintf(fout, '%s\n', sequence2);
fclose(fout);

RNApdistResult = tempname;
RunCmd('RNApdist', RNApdistPair, RNApdistResult)

delete(RNApdistPair);

fid = fopen(RNApdistResult, 'r');
P = fscanf(fid, '%f');
fclose(fid);

delete(RNApdistResult);

linenum = length(P);
RNApdist = (P(1:linenum));