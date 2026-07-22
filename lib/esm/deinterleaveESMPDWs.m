    function [labels,tracks] = deinterleaveESMPDWs(P,pdw,E)
        % Unsupervised weighted clustering followed by PRI-phase refinement.
        % Features combine RF, coarse DOA, PW and a local inter-pulse interval.
        m=numel(pdw); K=min(P.esmN,m); labels=[]; tracks=[];
        if m==0 || K==0, return; end
        toa=[pdw.toa].'; f=[pdw.freq].'; doa=[pdw.doa].'; pw=[pdw.pw].';
        lpri=localPRIfeature(toa,f,doa,pw,E);
        zpw=log(max(pw,eps));
        F=[robustScale(f)*0.65, robustScale(doa)*1.45, ...
            robustScale(zpw), robustScale(lpri)*1.25];
        labels=simpleKmeansLocal(F,K);

        % PRI-phase re-assignment gives the clusters pulse-train memory instead
        % of treating every PDW as an independent point.
        for it=1:3
            stats=trackStatsLocal(pdw,labels,K);
            new=labels;
            for i=1:m
                best=inf; bk=labels(i);
                for k=1:K
                    if stats(k).count==0, continue; end
                    c=1.3*abs(doa(i)-stats(k).doa)/max(3,stats(k).doaStd+2);
                    c=c+abs(log(max(pw(i),eps)/max(stats(k).pw,eps)))/0.28;
                    ff=stats(k).freqSamples;
                    c=c+0.65*min(abs(f(i)-ff))/max(0.16e6,0.18*(stats(k).freqSpan+0.2e6));
                    if stats(k).PRI>0
                        u=(toa(i)-stats(k).t0)/stats(k).PRI;
                        c=c+1.15*abs(u-round(u))/0.22;
                    end
                    if c<best, best=c; bk=k; end
                end
                new(i)=bk;
            end
            if all(new==labels), break; end
            labels=new;
        end

        stats=trackStatsLocal(pdw,labels,K);
        % Stable display order by median DOA.
        [~,ord]=sort([stats.doa]); remap=zeros(1,K);
        for k=1:K, remap(ord(k))=k; end
        labels=arrayfun(@(x)remap(x),labels);
        stats=trackStatsLocal(pdw,labels,K);
        tracks=stats;
        for k=1:K
            tracks(k).id=k;
            tracks(k).class=classifyESMTrack(tracks(k));
            ids=[pdw(labels==k).truthId];
            if isempty(ids), tracks(k).truthId=NaN; tracks(k).purity=0;
            else
                u=unique(ids); cnt=zeros(size(u));
                for q=1:numel(u), cnt(q)=sum(ids==u(q)); end
                [mx,ii]=max(cnt); tracks(k).truthId=u(ii); tracks(k).purity=mx/numel(ids);
            end
            tracks(k).range=P.c*tracks(k).t0;
        end
    end
