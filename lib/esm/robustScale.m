    function z = robustScale(x)
        x=x(:); med=median(x); sc=median(abs(x-med));
        if ~(sc>0), sc=(max(x)-min(x))/4; end
        if ~(sc>0), sc=1; end
        z=(x-med)/sc;
    end
