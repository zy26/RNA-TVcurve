function MultipleRNAP(input_path, seq_file, stru_file, output_path, path)
global globalpath
globalpath = path;

mkdir(fullfile(output_path, 'Related_results')) %% save detailed results

%%%%%%%%%  download file name: TV curve and similarity matrix
Simwaveletfile = fullfile(output_path, 'aaa-Simwavelet.txt');

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

RNAdistancePair = 'aaa-RNAdistancePair.txt';
fout1 = fopen(RNAdistancePair, 'w');

RNApdistPair = 'aaa-RNApdistPair.txt';
fout2 = fopen(RNApdistPair, 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%计算相似性

Sim_wavelet = zeros(N, N); %%%小波的多分辨率相似性：

for i = 1:N
    Sim_wavelet(i, i) = 1;
    numseq{i} = nummericalRNAother(seq{i}, stru{i});
    for j = i + 1:N
        %%%保存成 RNAdistance可以批处理的形式 aaa-RNApdistPair.txt
        %%%保存成 RNAdistance可以批处理的形式 aaa-RNAdistancePair.txt
     
        fprintf(fout1, '%s\n', stru{i});
        fprintf(fout1, '%s\n', stru{j});
        fprintf(fout2, '%s\n', seq{i});
        fprintf(fout2, '%s\n', seq{j});
     
        %%%%%%%%%%%%%%%%%%
        %%%小波的多分辨率相似性：
        Sim_wavelet(i, j) = waveletSimilarity(numseq{i}, numseq{j});
        Sim_wavelet(j, i) = Sim_wavelet(i, j);
    end
end
fclose(fout1);
fclose(fout2);


RNAdistanceResult = 'aaa-RNAdistanceResult.txt';
delete(RNAdistanceResult);
RunCmd('RNAdistance', RNAdistancePair, RNAdistanceResult);

%%%%%%%%%%%读取利用RNAdistance得到的距离文件

fid = fopen(RNAdistanceResult, 'r');
P = fscanf(fid, ['%s %f']);
fclose(fid);
linenum = length(P) / 3;
RNAdistance = squareform(P(linenum * 2 + 1:linenum * 3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RNApdistResult = 'aaa-RNApdistResult.txt';

delete(RNApdistResult);
RunCmd('RNApdist', RNApdistPair, RNApdistResult);

%
%      %%%%%%%%%%%读取利用RNAdistance得到的距离文件
%
fid = fopen(RNApdistResult, 'r');
P = fscanf(fid, ['%f']);
fclose(fid);
linenum = length(P);
RNApdist = squareform(P(1:linenum));

Label = Header;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  save distance matrix
% 1 - Sim_wavelet
seqName = cellstr(Header)';

T1 = array2table(1 - Sim_wavelet, 'RowNames', seqName)

writetable(T1, fullfile(output_path, 'Distance_matrix_RNATVcurve.csv'), 'WriteRowNames', true);

T2 = array2table(RNApdist, 'RowNames', seqName)

writetable(T2, fullfile(output_path, 'Distance_matrix_RNApdist.csv'), 'WriteRowNames', true);

T3 = array2table(RNAdistance, 'RowNames', seqName)

writetable(T3, fullfile(output_path, 'Distance_matrix_RNAdistance.csv'), 'WriteRowNames', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Build phylogenetic tree
Tree = seqneighjoin_average((1 - Sim_wavelet), 'average', Label);
phytreewrite(fullfile(output_path, 'tree_RNATVcurve.tree'), Tree)

if (N > 20)
    plot(Tree,'orient','top');
else
    plot(Tree);
end
set(gcf,'PaperPositionMode','auto')
print('-dsvg', fullfile(output_path, 'Tree_RNATVcurve'), '-r0')

try
    %%%%% Build phylogenetic tree
    Treedistance = seqneighjoin_average(RNAdistance, 'average', Label);
    phytreewrite(fullfile(output_path, 'tree_RNAdistance.tree'), Treedistance)
    if (N > 20)
        plot(Treedistance,'orient','top');
    else
        plot(Treedistance);
    end
    set(gcf,'PaperPositionMode','auto')
    print('-dsvg', fullfile(output_path, 'tree_RNAdistance'), '-r0')
catch
end

try
    %%%%% Build phylogenetic tree
    Treepdist = seqneighjoin_average(RNApdist, 'average', Label);
    phytreewrite(fullfile(output_path, 'tree_RNApdist.tree'), Treepdist)
    if (N > 20)
        plot(Treepdist,'orient','top');
    else
        plot(Treepdist);
    end
    set(gcf,'PaperPositionMode','auto')
    print('-dsvg', fullfile(output_path, 'tree_RNApdist'), '-r0')
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%  save distance matrix
% 1 - Sim_wavelet

save(Simwaveletfile, '-ascii', '-double', 'Sim_wavelet');
save(fullfile(output_path, 'aaa-RNAdistance.txt'), '-ascii', '-double', 'RNAdistance');
save(fullfile(output_path, 'aaa-RNApdist.txt'), '-ascii', '-double', 'RNApdist');

try
    movefile('*.txt', output_path)
catch
end

try
    delete('*.ps')
catch
end