function [header, numseq] = loadPairFile(id, input_path, seq_file, stru_file, output_path)

[header, seq] = fastaread(fullfile(input_path, seq_file));
N = length(seq);

fid = fopen(fullfile(input_path, stru_file));
stru = textscan(fid, '%s', 'Delimiter', '\n');
stru = stru{1};
fclose(fid);

numseq = cell(N, 1);

if iscell(seq)
    for i = 1:N
        mkdir(fullfile(output_path, sprintf('%03d', i)));
        numseq{i} = nummericalRNA(seq{i}, stru{i}, fullfile(output_path, sprintf('%03d', i), sprintf('RNA-TVcurve-%s', id)), header{i});
        RNAplot(seq{i}, stru{i}, fullfile(output_path, sprintf('%03d',i), sprintf('RNA-Structure-%s', id)));
        fout = fopen(fullfile(output_path, sprintf('%03d',i), sprintf('aaa-TVcurve-%s.txt', id)), 'w');
        fprintf(fout, '%s\n', strcat('>', header{i}));
        fprintf(fout, ' %3.0f\n', numseq{i});
        fclose(fout);
    end
else
    mkdir(fullfile(output_path));
    numseq = nummericalRNA(seq, stru{1}, fullfile(output_path, sprintf('RNA-TVcurve-%s', id)), header);
    RNAplot(seq, stru{1}, fullfile(output_path, sprintf('RNA-Structure-%s', id)));
    fout = fopen(fullfile(output_path, sprintf('aaa-TVcurve-%s.txt', id)), 'w');
    fprintf(fout, '%s\n', strcat('>', header));
    fprintf(fout, ' %3.0f\n', numseq);
    fclose(fout);
end
