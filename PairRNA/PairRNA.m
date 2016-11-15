%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of PairRNA.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PairRNA(input_path, seq_file_1, stru_file_1, seq_file_2, stru_file_2, output_path, global_path)
global globalpath
globalpath = global_path; % Maybe we can change global_path to matlab cd

%%%%%%%%%%%%%%%%%%%%%%  Read fileframe-1
[header1, numseq_1] = loadPairFile('1', input_path, seq_file_1, stru_file_1, output_path);

%%%%%%%%%%%%%%%%%%%%%%  Read fileframe-2

[header2, numseq_2] = loadPairFile('2', input_path, seq_file_2, stru_file_2, output_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RNA-TV的相似性度量

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%计算相似性
N = length(header1);

if (N ~= length(header2))
    ME = MException('MyComponent:LengthNotEqual', ...
        'Length Not Equal: %s %s', header1, header2);
    throw(ME);
end

Sim_wavelet = zeros(N, 1); %%%小波的多分辨率相似性：

if iscell(header1)
    for i = 1:N
        %%%%%%%%%%%%%%%%%%
        %%%小波的多分辨率相似性：
        Sim_wavelet(i) = waveletSimilarity(numseq_1{i}, numseq_2{i});
        savedouble(fullfile(output_path, sprintf('%03d', i), 'distance_wavelet.txt'), header1{i}, header2{i}, 1 - Sim_wavelet(i));
    end
else
    %%%%%%%%%%%%%%%%%%
    %%%小波的多分辨率相似性：
    Sim_wavelet = waveletSimilarity(numseq_1, numseq_2);
    savedouble(fullfile(output_path, 'distance_wavelet.txt'), header1, header2, 1 - Sim_wavelet);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of PairRNA_webfunction.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%