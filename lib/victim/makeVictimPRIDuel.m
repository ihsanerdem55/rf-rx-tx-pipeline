    function x=makeVictimPRIDuel(k,replay,ref,guard,Nfast,Fs,skinAmp, ...
            ghostEff,fdPulse,phiPulse,victimHop,coverDelays,noisePow,rs,P)
        L=numel(ref); tp=(0:L-1)/Fs; x=complex(zeros(Nfast,1));
        x=addPulseAtDelay(x,skinAmp*ref(:),0,guard,Fs);
        for m=1:size(ghostEff,1)
            d=ghostEff(m,k);
            mismatch=exp(-1i*2*pi*victimHop(k)*tp);
            dop=exp(1i*(phiPulse(k)+2*pi*fdPulse(k)*tp));
            sig=1.15*replay(:).*mismatch(:).*dop(:);
            x=addPulseAtDelay(x,sig,d,guard,Fs);
        end
        if P.cover
            for q=1:numel(coverDelays)
                jitter=(rand(rs)-0.5)*0.20*P.PW;
                ph=2*pi*rand(rs);
                sig=P.coverAmp*replay(:).*exp(1i*ph);
                x=addPulseAtDelay(x,sig,max(0,coverDelays(q)+jitter),guard,Fs);
            end
        end
        if ~strcmp(P.jammode,'Off')
            jam=victimJammer(Nfast,Fs,P,victimHop(k),rs);
            jamPow=(skinAmp^2/max(L,1))*10^(P.jsdb/10);
            jam=jam/sqrt(mean(abs(jam).^2)+eps)*sqrt(jamPow);
            x=x+jam(:);
        end
        x=x+sqrt(noisePow/2)*(randn(rs,Nfast,1)+1i*randn(rs,Nfast,1));
    end
