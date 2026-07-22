    function y = firZeroPhase(x, hf)
        % A Type-I FIR has exactly linear phase, i.e. a constant group delay of
        % M = (Ntaps-1)/2 samples. Convolve, then drop exactly that delay: the
        % result is truly zero-phase with NO magnitude squaring.
        x = x(:).'; N = numel(x);
        M = (numel(hf)-1)/2;
        yfull = conv(x, hf);
        y = yfull(M+1 : M+N);
    end
