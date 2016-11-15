function ScoreAnalysis(Score, Header, Sequence, Name, LL)
%% Score:N*3;
[~, order] = sort(Score);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sequence, Structure? Header???????????
order = order + 1; %%%%%%%%%%%%%%%%%%%
orderHeader = Header(order);
%%%%% Look the first LL maximal difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:LL
    newrnaps = strcat(Name, num2str(i), orderHeader{i});
    SingleRNAfold(Sequence{(order(i))}, newrnaps);
end