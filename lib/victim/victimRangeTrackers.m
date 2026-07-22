    function [peakTrack,leadTrack]=victimRangeTrackers(rc,rng,P,res,maxDelay)
        Np=size(rc,2); peakTrack=nan(1,Np); leadTrack=nan(1,Np);
        lo=P.Rng-res; hi=P.Rng+P.c*maxDelay/2+4*res;
        gate=find(rng>=lo & rng<=hi); if isempty(gate), gate=1:numel(rng); end
        for k=1:Np
            a=abs(rc(:,k)); [~,ip]=max(a(gate)); peakTrack(k)=rng(gate(ip));
            med=median(a); sig=1.4826*median(abs(a-med))+eps;
            thr=med+4*sig; q=gate(a(gate)>thr);
            if isempty(q), leadTrack(k)=peakTrack(k); else, leadTrack(k)=rng(q(1)); end
        end
    end
