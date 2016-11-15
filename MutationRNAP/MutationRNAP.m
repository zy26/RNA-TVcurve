function MutationRNAP(input_path, seq_file, output_path, global_path)
global globalpath
globalpath = global_path; % Maybe we can change global_path to matlab cd

mkdir(fullfile(output_path, 'Related_results')) %% save detailed results
%%%%%%%%%%%%%%%%%% save jpg
deleteriousProfile_fileName = fullfile(output_path, 'deleteriousProfile');
deleteriousHist_fileName = fullfile(output_path, 'deleteriousHist');

%%%%%%%%%%%%%%%%%%%%%%
seq = ReadSequenceFile(fullfile(input_path, seq_file));

Seq = dna2rna(seq);
Pos = length(Seq);

[Header, Sequence, Structure, Energe] = fMutationAllRNAFoldFile(Seq, 'Original');

%%%%% the first in Sequence,Structure,Energe is the original Seq
MutationNum = length(Header) - 1;

numOrigSeq = nummericalRNA(Sequence{1}, Structure{1});

TVcurve = cell(1, MutationNum);

TVcurve{1} = numOrigSeq;

SimTV_wavelet = zeros(MutationNum, 1);
SimRNAdistance = zeros(MutationNum, 1);
SimRNApdist = zeros(MutationNum, 1);

numEnerge = zeros(1, MutationNum);
numMutationSeq = cell(1, MutationNum);
deleteriousSim = zeros(1, fix(MutationNum / 3));
deleteriousSimRNAdistance = zeros(1, fix(MutationNum / 3));
deleteriousSimRNApdist = zeros(1, fix(MutationNum / 3));
deleteriousEnergy = zeros(1, fix(MutationNum / 3));
deleteriousRNAdistanceHeader = cell(1, fix(MutationNum / 3));
deleteriousHeader = cell(1, fix(MutationNum / 3));
deleteriousRNApdistHeader = cell(1, fix(MutationNum / 3));

for i = 1:MutationNum
    numEnerge(i) = str2double(Energe{i + 1});
    numMutationSeq{i} = nummericalRNA(Sequence{i + 1}, Structure{i + 1});
    TVcurve{i + 1} = numMutationSeq{i};
    SimTV_wavelet(i) = waveletSimilarity(numOrigSeq, numMutationSeq{i});
 
    SimRNAdistance(i) = RNADistance(Structure{1}, Structure{i + 1});
    SimRNApdist(i) = RNApdist(Sequence{1}, Sequence{i + 1});
 
    if rem(i, 3) == 0
        t = fix(i / 3);
        %% mean deterious profile
        [deleteriousSim(t)] = mean(SimTV_wavelet(i - 2:i));
        [deleteriousSimRNAdistance(t)] = mean(SimRNAdistance(i - 2:i));
        [deleteriousSimRNApdist(t)] = mean(SimRNApdist(i - 2:i));
     
        %% similar to RDMAS deleterious pofile
        [deleteriousSim(t), pp] = max(SimTV_wavelet(i - 2:i));
        [deleteriousSimRNAdistance(t), ppRNAdistance] = max(SimRNAdistance(i - 2:i));
        [deleteriousSimRNApdist(t), ppRNApdist] = max(SimRNApdist(i - 2:i));
        deleteriousEnergy(t) = numEnerge(i - 3 + pp);
        deleteriousRNAdistanceHeader{t} = Header{1 + i - 3 + ppRNAdistance};
        deleteriousHeader{t} = Header{1 + i - 3 + pp};
        deleteriousRNApdistHeader{t} = Header{1 + i - 3 + ppRNApdist};
    end
end

SimTV_wavelet = 1 - SimTV_wavelet; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, order] = sort(SimTV_wavelet, 'descend');
disp(['the maximal mutation of TV-curve is', Header{order(1) + 1}]);

%%%%%%%%%%%%%%%%%%%%%%%
[~, orderRNAdistance] = sort(SimRNAdistance, 'descend');
disp(['the maximal mutation of RNAdistance is  ', Header{orderRNAdistance(1) + 1}]);

%%%%%%%%%%%%%%%%%%%%%%%%
[~, orderRNApdist] = sort(SimRNApdist, 'descend');
disp(['the maximal mutation of RNApdist', Header{orderRNApdist(1) + 1}]);

%     TvOrder_Head= Header{order + 1}
%     PdistOrder_Head= Header{orderRNApdist + 1}
%     DistanceOrder_Head=Header{orderRNAdistance + 1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Save Detail file-1:mutated sequence and structures,TV-curve

if(exist(fullfile(output_path, 'Related_results', 'MinimumFreeEnergy_all.txt'), 'file'))
    delete(fullfile(output_path, 'Related_results', 'MinimumFreeEnergy_all.txt'));
end
if(exist(fullfile(output_path, 'Related_results', 'Sequence_all.txt'), 'file'))
    delete(fullfile(output_path, 'Related_results', 'Sequence_all.txt'));
end
if(exist(fullfile(output_path, 'Related_results', 'Structure_all.txt'), 'file'))
    delete(fullfile(output_path, 'Related_results', 'Structure_all.txt'));
end
fastawrite(fullfile(output_path, 'Related_results', 'MinimumFreeEnergy_all.txt'), Header, Energe)
fastawrite(fullfile(output_path, 'Related_results', 'Sequence_all.txt'), Header, Sequence)
fastawrite(fullfile(output_path, 'Related_results', 'Structure_all.txt'), Header, Structure)

fout2 = fopen(fullfile(output_path, 'Related_results', 'TVcurve_all.txt'), 'w');
for i = 1:MutationNum + 1
    fprintf(fout2, '%s\n', strcat('>', Header{i}));
    fprintf(fout2, '%5.1f\n', TVcurve{i});
end
fclose(fout2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Save Detail file-2:Rank_MutatedStructure_for3methods

oHeader = cellstr(Header);

Rank = (1:length(order))';
results_RNATVcurve = (oHeader(order + 1))';
results_RNApdist = (oHeader(orderRNApdist + 1))';
results_RNAdistance = (oHeader(orderRNAdistance + 1))';

ALL = table(Rank, results_RNATVcurve, results_RNApdist, results_RNAdistance);
writetable(ALL, fullfile(output_path, 'Related_results', 'Rank_MutatedStructure_for3methods.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SingleRNAfold(Sequence{1}, fullfile(output_path, 'original.txt'));
SingleRNAfold(Sequence{order(1) + 1}, fullfile(output_path, strcat('RNA-TVcurve_', Header{order(1) + 1}, '.txt')));
SingleRNAfold(Sequence{orderRNAdistance(1) + 1}, fullfile(output_path, strcat('RNAdistance_', Header{orderRNAdistance(1) + 1}, '.txt')));
SingleRNAfold(Sequence{orderRNApdist(1) + 1}, fullfile(output_path, strcat('RNApdist_', Header{orderRNApdist(1) + 1}, '.txt')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the maximum  mutation on each position
disp('correlation analysis between 1-deleteriousSim and deleteriousSimRNAdistance');
corr(1 - deleteriousSim',deleteriousSimRNAdistance')
disp('correlation analysis between 1-deleteriousSim and deleteriousSimRNApdist');
corr(1 - deleteriousSim',deleteriousSimRNApdist')

clf
x = 1:Pos;
deleteriousSim = (1 - deleteriousSim - mean(1 - deleteriousSim)) / std(1 - deleteriousSim);
deleteriousSimRNAdistance = (deleteriousSimRNAdistance - mean(deleteriousSimRNAdistance)) / std(deleteriousSimRNAdistance);
deleteriousSimRNApdist = (deleteriousSimRNApdist - mean(deleteriousSimRNApdist)) / std(deleteriousSimRNApdist);
plot(x, deleteriousSim, 'r', x, deleteriousSimRNAdistance, 'b', x, deleteriousSimRNApdist, 'g');
h = legend('RNA-TVcurve', 'RNAdistance', 'RNApdist', 1);
set(h);
print('-dsvg', deleteriousProfile_fileName)

clf

subplot(131), hist(deleteriousSim); title('RNA-TVcurve')
subplot(132), hist(deleteriousSimRNAdistance); title('RNAdistance')
subplot(133), hist(deleteriousSimRNApdist); title('RNApdist')
print('-dsvg', deleteriousHist_fileName)

%%% original
SingleRNAfold(Sequence{1}, fullfile(output_path, 'Original.txt'));

delete('1_dp.ps');
delete('2_dp.ps');
%movefile('0.ps', 'Original secondary structure.ps');