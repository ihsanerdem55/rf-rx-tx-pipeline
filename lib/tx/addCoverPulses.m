    function [y,delays] = addCoverPulses(y,base,P,seed)
        % A bank of lower-power DRFM pulse trains distributed across the range
        % gate. Random fixed phases prevent them from adding as one ideal target.
        delays=[];
        if ~P.cover || isempty(base), return, end
        N=max(1,round(P.coverN));
        span=min(0.8*P.PRI,max([4*P.PW,P.drdelay+2*P.PW,2*P.PW]));
        delays=linspace(0.15*P.PW,span,N);
        dmax=round(max(delays)*P.Fs);
        wasReal=isreal(y);
        if wasReal
            zy=analytic(y); zb=analytic(base);
        else
            zy=y; zb=base;
        end
        if numel(zy)<numel(base)+dmax, zy(end+1:numel(base)+dmax)=0; end
        rs=RandStream('mt19937ar','Seed',seed);
        for q=1:N
            d=round(delays(q)*P.Fs);
            ph=2*pi*rand(rs);
            seg=d+(1:numel(zb));
            zy(seg)=zy(seg)+P.coverAmp*zb.*exp(1i*ph);
        end
        if wasReal, y=real(zy); else, y=zy; end
    end
