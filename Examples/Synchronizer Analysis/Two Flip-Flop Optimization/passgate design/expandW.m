function resW = expandW(result)
%expandW takes the result passed into this function and expands the result
%to the appropriate size for the synchronizer in question
%   This expandW function is for a two flip-flop passgate synchronizer

bufN = 1;
bufP = 7;
latchN = 2:6;
latchP = 8:length(result);

resW = [result(bufN);
        result(latchN);
        result(latchN);
        result(latchN);
        result(latchN);
        result(bufP);
        result(latchP);
        result(latchP);
        result(latchP);
        result(latchP)];

end
