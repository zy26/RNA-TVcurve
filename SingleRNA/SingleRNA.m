%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of SingleRNA.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SingleRNA(input_path, seq_file, stru_file, output_path, global_path)
global globalpath
globalpath = global_path; % Maybe we can change global_path to matlab cd

[header, seq] = fastaread(fullfile(input_path, seq_file));

fid = fopen(fullfile(input_path, stru_file));
stru = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
stru = stru{1};

if iscell(seq)
    for i = 1:length(seq)
        mkdir(fullfile(output_path, sprintf('%03d',i)));
        numseq = nummericalRNA(seq{i}, stru{i}, fullfile(output_path, sprintf('%03d',i), 'aaa-TV-Curve'), header{i});
        waveletTransformPlot(numseq, fullfile(output_path, sprintf('%03d',i), '1 TV-curve and approximation singnal of wavelet decomposition'), fullfile(output_path, sprintf('%03d',i), '2 Detailed signals of wavelet-decomposition'));
        RNAplot(seq{i}, stru{i}, fullfile(output_path, sprintf('%03d',i), '0 RNA structure'));
        SaveStructure(header{i}, seq{i}, stru{i}, fullfile(output_path, sprintf('%03d',i), 'Sequence structure'));
    end
else
    numseq = nummericalRNA(seq, stru{1}, fullfile(output_path, 'aaa-TV-Curve'), header);
    waveletTransformPlot(numseq, fullfile(output_path, '1 TV-curve and approximation singnal of wavelet decomposition'), fullfile(output_path, '2 Detailed signals of wavelet-decomposition'));
    RNAplot(seq, stru{1}, fullfile(output_path, '0 RNA structure'));
    SaveStructure(header, seq, stru{1}, fullfile(output_path, 'Sequence structure'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of SingleRNA_webfunction.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%