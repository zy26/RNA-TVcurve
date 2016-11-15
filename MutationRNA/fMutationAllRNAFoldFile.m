
%>P5abc subdomain of the group I intron ribozyme in Tetrahymena thermophila

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % >Lep RNA
% % >The spliced leader RNA gene of Leptomonas collosoma
% % AACTAAAACAATTTTTGAAGAACAGTTTCTGTACTTCATTGGTATGTAGAGACTTC

% Seqname='Lep'
% Lep='AACTAAAACAATTTTTGAAGAACAGTTTCTGTACTTCATTGGTATGTAGAGACTTC';
% Lep=wrev(Lep)
% Seq=Lep;

function [Header, Sequence, Structure, Energe] = fMutationAllRNAFoldFile(Seq, Seqname)

%%%%%%%%%%%%%%%%  for the  the mutation of  all the Seq

Seq = dna2rna(Seq);

MutationRNAfoldSeq = tempname;

fout = fopen(MutationRNAfoldSeq, 'w');

len = length(Seq);
%%%%%%%%%%%%%%%%%%%%  Using RNAfold.exe £¨MutationRNAfoldSeq.txt£© to fold all single point mutation
%%%%%%%%%%%%%%%%%%%%  RNA sequence and save the results in the txt file of
%%%%%%%%%%%%%%%%%%%%  MutationRNAfoldSeqStru.txt
%%%%%%%%%%%%%%%%%%%% the first of  MutationRNAfoldSeq.txt and
%%%%%%%%%%%%%%%%%%%% MutationRNAfoldSeqStru.txt is the original sequence
mutation = ['A', 'U', 'G', 'C'];
%%%% save the original sequence in the file   MutationRNAfoldSeq.txt
fprintf(fout, '%s', Seq);
fprintf(fout, '\r\n');
MutationseqCell = cell(1, len * 4 + 1);
MutationHeaderCell = cell(1, len * 4 + 1);
MutationseqCell{1} = Seq;
MutationHeaderCell{1} = Seqname;
t = 2;
for i = 1:len
    for j = 1:4
        MutationSeq = Seq;
        if ~ strcmp(MutationSeq(i), mutation(j))
            MutationSeq(i) = mutation(j);
            fprintf(fout, '%s', MutationSeq);
            fprintf(fout, '\r\n');
            MutationseqCell{t} = MutationSeq;
            MutationHeaderCell{t} = strcat(Seq(i), num2str(i), mutation(j));
            t = t + 1;
        end
    end
end

MutationseqCell(t:len * 4 + 1) = [];
MutationHeaderCell(t:len * 4 + 1) = [];

fclose(fout);

MutationRNAfoldSeqStru = tempname;
RunCmd('RNAfold', MutationRNAfoldSeq, MutationRNAfoldSeqStru)
delete(MutationRNAfoldSeq);

Header = MutationHeaderCell;
Sequence = MutationseqCell;

Num = length(Header);
Structure = cell(1, Num);
Energe = cell(1, Num);
fin = fopen(MutationRNAfoldSeqStru);
tt = 0;
nn = 0;
tline = fgetl(fin);
while ischar(tline)
    tt = tt + 1;
    if rem(tt, 2) == 0
        nn = tt / 2;
        len = length(Sequence{nn});
        len2 = length(tline);
        Structure{nn} = tline(1:len);
        pp = tline(len + 1:len2);
        pp = strrep(pp, '(', '');
        pp = strrep(pp, ')', '');
        Energe{nn} = pp;
    end
    tline = fgetl(fin);
end
Structure(nn + 1:Num) = [];
Energe(nn + 1:Num) = [];

fclose(fin);
delete(MutationRNAfoldSeqStru);
