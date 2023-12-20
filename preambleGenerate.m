function [preamble] = preambleGgenerate(length)
% preambleGenerate() 
% input : length: a scaler value, desired length of preamble.
% output: preamble: preamble bits
    preamble = zeros(length, 1);

    lfsr = ones(8);

    for ii = 1:length
        preamble(ii) = lfsr(8);
        lfsrInter = bitxor(lfsr(8),lfsr(6));
        lfsrInter = bitxor(lfsrInter,lfsr(5));
        lfsrInter = bitxor(lfsrInter, lfsr(4));
        lfsr = circshift(lfsr,1);
        lfsr(1) = lfsrInter;
    end
end
