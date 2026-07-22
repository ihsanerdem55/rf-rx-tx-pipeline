    function x=addPulseAtDelay(x,sig,dly,guard,Fs)
        sig=sig(:); L=numel(sig); d=max(0,round(dly*Fs));
        idx=guard+d+(1:L); good=idx<=numel(x);
        if any(good), x(idx(good))=x(idx(good))+sig(good); end
    end
