    function pdw = extractESMPDWs(P,X,t,truth,noiseStd,seed)
        % Educational PDW front end: pulse TOA/PW/frequency are measured with
        % finite errors, while a conventional three-element beam scan supplies
        % a deliberately coarse DOA. MUSIC later replaces this rough angle.
        d0=struct('toa',0,'pw',0,'freq',0,'doa',0,'amp',0,'range',0, ...
            'truthId',0,'snr',0);
        pdw=repmat(d0,0,1); rs=RandStream('mt19937ar','Seed',seed);
        pos=P.esmPos(:); grid=-90:1:90;
        for k=1:numel(truth)
            q=truth(k);
            snr=20*log10(max(q.amp,eps)/max(sqrt(2)*noiseStd,eps));
            miss=P.pdwMiss;
            for j=1:numel(truth)
                if j==k, continue; end
                z=truth(j);
                overlap=(q.toa<z.toa+z.pw) && (z.toa<q.toa+q.pw);
                sameChan=abs(q.freq-z.freq)<max(0.20e6,0.6/max(q.pw,z.pw));
                if overlap && sameChan && z.amp>q.amp, miss=miss+0.22; end
            end
            miss=miss+0.10*max(0,(6-snr)/12);
            if rand(rs)<min(0.75,miss), continue; end

            toa=q.toa+randn(rs)*max(1/P.Fs,0.012*q.pw);
            pw=max(2/P.Fs,q.pw*(1+0.035*randn(rs)));
            f=q.freq+randn(rs)*max(4e3,0.015/max(q.pw,eps));
            c0=toa+0.5*pw;
            half=max(3,round(min(0.35*pw,0.8e-6)*P.Fs/2));
            ic=round(c0*P.Fs)+1; idx=max(1,ic-half):min(numel(t),ic+half);
            if isempty(idx), continue; end
            v=mean(X(:,idx).*exp(-1i*2*pi*f*t(idx)),2);
            rf=P.Flo+f; lam=P.c/max(rf,eps);
            sc=zeros(size(grid));
            for a=1:numel(grid)
                sv=exp(1i*2*pi*pos*sind(grid(a))/lam);
                sc(a)=abs(sv'*v)^2/(real(sv'*sv)*real(v'*v)+eps);
            end
            [~,im]=max(sc); doa=grid(im)+0.25*randn(rs);
            d=d0; d.toa=max(0,toa); d.pw=pw; d.freq=f; d.doa=doa;
            d.amp=norm(v); d.range=P.c*d.toa; d.truthId=q.id; d.snr=snr;
            pdw(end+1)=d; %#ok<AGROW>
        end
        if ~isempty(pdw), [~,ix]=sort([pdw.toa]); pdw=pdw(ix); end
    end
