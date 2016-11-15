function MultipleRNA(input_path, seq_file, stru_file, output_path)

%%%%%%%%%%%%%%%%%%%%%%
[Header, seq] = fastaread(fullfile(input_path, seq_file));
N = length(seq);

fid = fopen(fullfile(input_path, stru_file));
stru = textscan(fid, '%s', 'Delimiter', '\n');
stru = stru{1};
fclose(fid);

fout = fopen(fullfile(output_path, 'aaa-TVcurve.txt'), 'w');
fout1 = fopen(fullfile(output_path, 'aaa-Sequence.txt'), 'w');
fout2 = fopen(fullfile(output_path, 'aaa-Structure.txt'), 'w');
numseq = cell(N, 1);
for i = 1:N
    numseq{i} = nummericalRNAother(seq{i}, stru{i});
    fprintf(fout, '%s\n', strcat('>', Header{i}));
    fprintf(fout, ' %3.0f', numseq{i});
    fprintf(fout, '\n');
    fprintf(fout2, '%s\n', strcat('>', Header{i}));
    fprintf(fout2, '%s\n', stru{i});
    fprintf(fout1, '%s\n', strcat('>', Header{i}));
    fprintf(fout1, '%s\n', seq{i});
end
fclose(fout);
fclose(fout1);
fclose(fout2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%计算相似性
Sim_wavelet = zeros(N, N); %%%小波的多分辨率相似性：

for i = 1:N
    Sim_wavelet(i, i) = 1;
    for j = i + 1:N
        %%%%%%%%%%%%%%%%%%
        %%%小波的多分辨率相似性：
        Sim_wavelet(i, j) = waveletSimilarity(numseq{i}, numseq{j});
        Sim_wavelet(j, i) = Sim_wavelet(i, j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%  save distance matrix
% 1 - Sim_wavelet
seqName = cellstr(Header)';
T1 = array2table(1 - Sim_wavelet, 'RowNames', seqName);

writetable(T1, fullfile(output_path, 'Distance_matrix_RNATVcurve.csv'), 'WriteRowNames', true);

Label = Header;
%%%% Build phylogenetic tree
Tree = seqneighjoin_average((1 - Sim_wavelet), 'average', Label);
phytreewrite(fullfile(output_path, 'tree_RNATVcurve.tree'), Tree)
if (N > 20)
    plot(Tree,'orient','top');
else
    plot(Tree);
end
set(gcf,'PaperPositionMode','auto')
print('-dsvg', fullfile(output_path, 'Tree_RNATVcurve'), '-r0')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    movefile('*.txt', output_path)
catch
end

try
    delete('*.ps')
catch
end