    function rc=addVictimClutter(rc,rangeAxis,tSlow,skinAmp,P)
        Nr=size(rc,1); rs=RandStream('mt19937ar','Seed',P.seed+9090);
        g=(randn(rs,Nr,1)+1i*randn(rs,Nr,1))/sqrt(2);
        ker=ones(9,1)/9; g=conv(g,ker,'same');
        roll=1./sqrt(1+max(rangeAxis(:)-min(rangeAxis),0)/max(P.Rng,1));
        g=g.*roll; g=g/sqrt(mean(abs(g).^2)+eps);
        amp=skinAmp*10^(P.clutterDB/20);
        lambda=P.c/max(P.Ftrf,eps); fdStd=2*P.clutterV/lambda;
        fd=fdStd*randn(rs,Nr,1);
        C=amp*g.*exp(1i*2*pi*fd*tSlow);
        rc=rc+C;
    end
