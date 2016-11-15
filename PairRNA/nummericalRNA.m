%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of nummericalRNA.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function numseq = nummericalRNA(seq, stru, TvcurvefileallName, header)

seqstru = characterseq(seq, stru);

% 8 charchers sequerences: AUGC ,BVHD
len = length(seq);

x = 0:1:3 * len;
numseq = zeros(1, 3 * len);

for i = 1:len
    % 3 * i - 2
    if i == 1
        if seqstru(i) == 'A' || seqstru(i) == 'H' || seqstru(i) == 'D' || seqstru(i) == 'U'
            numseq(3 * i - 2) = 1;
        else
            numseq(3 * i - 2) = - 1;
        end
    else
        if seqstru(i) == 'A' || seqstru(i) == 'H' || seqstru(i) == 'D' || seqstru(i) == 'U'
            numseq(3 * i - 2) = numseq(3 * i - 3) + 1;
        else
            numseq(3 * i - 2) = numseq(3 * i - 3) - 1;
        end
    end
    % 3 * i - 1
    if seqstru(i) == 'A' || seqstru(i) == 'H' || seqstru(i) == 'V' || seqstru(i) == 'C'
        numseq(3 * i - 1) = numseq(3 * i - 2) + 1;
    else
        numseq(3 * i - 1) = numseq(3 * i - 2) - 1;
    end
    % 3 * i
    if seqstru(i) == 'A' || seqstru(i) == 'B' || seqstru(i) == 'V' || seqstru(i) == 'U'
        numseq(3 * i) = numseq(3 * i - 1) + 1;
    else
        numseq(3 * i) = numseq(3 * i - 1) - 1;
    end
 
end

numseq = [0, numseq];

clf
plot(x, numseq);
title(header,'Interpreter','none');
print('-dsvg', TvcurvefileallName)

function seqstru = characterseq(seq, stru)

len = length(seq);
length(stru);
seqstru = zeros(1, len);
for i = 1:len
    if stru(i) == '(' || stru(i) == ')'
        seqstru(i) = repla(seq(i));
    else
        seqstru(i) = seq(i);
    end
end

function s = repla(s1)
if s1 == 'A'
    s = 'B';
elseif s1 == 'U'
    s = 'V';
elseif s1 == 'G'
    s = 'H';
elseif s1 == 'C'
    s = 'D';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of nummericalRNA.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
