function SingleRNAfold(Sequence, outputFile)

singleRNAfoldSeq = tempname;

WriteSequence(singleRNAfoldSeq, Sequence);

RunCmd('RNAfold', singleRNAfoldSeq, outputFile);
RunCmd('RNAplot -o svg', outputFile, 'NUL');

delete(singleRNAfoldSeq);

%% renmae the output ps figure 'rna.ps' to '<File>.ps'
[pathstr,name,~] = fileparts(outputFile);

movefile('rna.svg', fullfile(pathstr,[name '.svg']))