    function [X,t,truth,noiseStd] = simulateESMScene(P,E,seed)
        % Complex analytic array capture. Using complex samples internally keeps
        % co-channel phase information intact even when the main GUI is set to
        % a real IF presentation.
        maxPRI=max([E.PRI]);
        Tobs=max(P.Np*P.PRI,11*maxPRI);
        Tobs=min(max(Tobs,0.9e-3),4.0e-3);
        maxDelay=max([E.range])/P.c;
        N=ceil((Tobs+maxDelay+max([E.PW])+8/P.Fs)*P.Fs);
        t=(0:N-1)/P.Fs; X=complex(zeros(3,N)); active=false(1,N);
        pos=P.esmPos(:);
        tt=struct('id',0,'toa',0,'emitTime',0,'pw',0,'freq',0,'rf',0, ...
            'az',0,'range',0,'amp',0,'classTruth','','pulseIndex',0,'pri',0);
        truth=repmat(tt,0,1);
        rs=RandStream('mt19937ar','Seed',seed);

        for e=1:numel(E)
            [emit,pris,hops]=esmEmitterSchedule(E(e),Tobs,rs);
            phi0=2*pi*rand(rs);
            for q=1:numel(emit)
                toa=emit(q)+E(e).range/P.c;
                idx=find(t>=toa & t<toa+E(e).PW);
                if isempty(idx), continue; end
                active(idx)=true;
                f=E(e).fc+hops(q)+E(e).vrad*E(e).rf/P.c;
                rf=P.Flo+f; lambda=P.c/max(rf,eps);
                tl=t(idx)-(toa+E(e).PW/2);
                intra=zeros(size(tl));
                if strcmp(E(e).modtype,'LFM chirp') && E(e).bchirp>0
                    intra=pi*(E(e).bchirp/E(e).PW)*tl.^2;
                end
                rangeGain=sqrt(max(P.Rng,1)/max(E(e).range,1));
                a=max(0.28,min(1.8,E(e).amp*rangeGain));
                b=a*exp(1i*(2*pi*f*t(idx)+intra+phi0));
                steer=exp(1i*2*pi*pos*sind(E(e).az)/lambda);
                X(:,idx)=X(:,idx)+steer*b;
                p=tt; p.id=e; p.toa=toa; p.emitTime=emit(q); p.pw=E(e).PW;
                p.freq=f; p.rf=rf; p.az=E(e).az; p.range=E(e).range;
                p.amp=a; p.classTruth=E(e).classTruth; p.pulseIndex=q;
                if q<=numel(pris), p.pri=pris(q); else, p.pri=E(e).PRI; end
                truth(end+1)=p; %#ok<AGROW>
            end
        end

        if any(active), Ps=mean(abs(X(1,active)).^2); else, Ps=1; end
        noiseStd=sqrt(Ps/max(10^(P.SNR/10),eps)/2);
        for ch=1:3
            X(ch,:)=X(ch,:)+noiseStd*(randn(rs,1,N)+1i*randn(rs,1,N));
        end
        if ~isempty(truth), [~,ix]=sort([truth.toa]); truth=truth(ix); end
    end
