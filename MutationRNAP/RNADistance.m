function [RNAdistance] = RNADistance(Structure1, Structure2)

RNAdistancePair = tempname;
fout = fopen(RNAdistancePair, 'w');
fprintf(fout, '%s\n', Structure1);
fprintf(fout, '%s\n', Structure2);
fclose(fout);

RNAdistanceResult = tempname;
RunCmd('RNAdistance', RNAdistancePair, RNAdistanceResult)

delete(RNAdistancePair);

fin = fopen(RNAdistanceResult, 'r');
P = fscanf(fin, ['%s %f']);
fclose(fin);

delete(RNAdistanceResult);

linenum = length(P) / 3;
RNAdistance = (P(linenum * 2 + 1:linenum * 3));