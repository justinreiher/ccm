function resW = expandW(result)
%expandW takes the result passed into this function and expands the result
%to the appropriate size for the synchronizer in question
%   This expandW function is for a two flip-flop jamb latch synchronizer

bufInN = 1;
bufOutN = 2;
bufInP = 8;
bufOutP = 9;
latchN = 3:7;
latchP = 10:length(result);

resW = [result(bufInN);
        result(bufOutN);
        result(latchN);
        result(latchN);
        result(latchN);
        result(latchN);
        result(bufInP);
        result(bufOutP);
        result(latchP);
        result(latchP);
        result(latchP);
        result(latchP)];

end

