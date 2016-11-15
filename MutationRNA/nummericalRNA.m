%%   http://rfam.janelia.org/browse.html

%%%%%
% (1,1),(1,1),(1,1) A\\
% (1,-1),(1,-1),(1,-1) A'\\
% (1,1),(1,1),(1,-1)U\\
% (1,-1),(1,-1),(1,1) U'\\
% (1,-1),(1,1),(1,1) G\\
% (1,1),(1,-1),(1,-1) G'\\
% (1,1),(1,-1),(1,1) C\\
% (1,-1),(1,1),(1,-1) C'
function numseq = nummericalRNA(seq, stru)
if nargin == 1, seqstru = seq;
end
if nargin == 2, seqstru = characterseq(seq, stru);
end
%8个字符的序列 A UGC ,B V H D
len = length(seq);

x = 0:1:3 * len;
numseq = zeros(1, 3 * len);

for i = 1:len
    %%%% 3*i-2
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
    %%%%% 3*i-1
    if seqstru(i) == 'A' || seqstru(i) == 'H' || seqstru(i) == 'V' || seqstru(i) == 'C'
        numseq(3 * i - 1) = numseq(3 * i - 2) + 1;
    else
        numseq(3 * i - 1) = numseq(3 * i - 2) - 1;
    end
    %%%%% 3*i
    if seqstru(i) == 'A' | seqstru(i) == 'B' | seqstru(i) == 'V' | seqstru(i) == 'U'
        numseq(3 * i) = numseq(3 * i - 1) + 1;
    else
        numseq(3 * i) = numseq(3 * i - 1) - 1;
    end
 
end

numseq = [0, numseq];

plot(x, numseq);

% seq = 'ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA';
% stru = '..(((((...((....))...))))).....';

function seqstru = characterseq(seq, stru)
%%%%%% 将二级结构序列转化为 8个字符的序列 A UGC ,B V H D

len = length(seq);
length(stru);
seqstru = zeros(len);
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