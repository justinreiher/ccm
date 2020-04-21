function resW = expandW(result)
%expandW takes the result passed into this function and expands the result
%to the appropriate size for the synchronizer in question
%   This expandW function is for a two flip-flop strongARM synchronizer

bufN1 = 1;
bufN2 = 2;
bufP1 = 13;
bufP2 = 14;
FFN = 3:12;
FFP = 15:length(result);

resW = [result(bufN1);
        result(bufN2);
        result(FFN);
        result(FFN);
        result(bufP1);
        result(bufP2);
        result(FFP);
        result(FFP);];

end

