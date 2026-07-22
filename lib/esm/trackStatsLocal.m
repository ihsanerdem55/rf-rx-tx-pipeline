    function st = trackStatsLocal(pdw,labels,K)
        T=struct('id',0,'count',0,'t0',0,'PRI',0,'priStd',0,'PW',0,'pw',0, ...
            'freq',0,'freqMed',0,'freqSpan',0,'freqLevels',0,'freqSamples',0, ...
            'doa',0,'doaStd',0,'range',0,'class','','truthId',NaN,'purity',0, ...
            'priAgile',false,'freqAgile',false,'agility','Fixed');
        st=repmat(T,1,K);
        for k=1:K
            idx=find(labels==k); st(k).id=k; st(k).count=numel(idx);
            if isempty(idx), st(k).freqSamples=0; continue; end
            tt=sort([pdw(idx).toa]); pp=[pdw(idx).pw]; ff=[pdw(idx).freq]; aa=[pdw(idx).doa];
            st(k).t0=min(tt); st(k).pw=median(pp); st(k).PW=st(k).pw;
            st(k).freq=median(ff); st(k).freqMed=st(k).freq;
            st(k).freqSpan=max(ff)-min(ff); st(k).freqSamples=ff;
            st(k).freqLevels=countFreqLevels(ff,0.12e6);
            st(k).doa=median(aa); st(k).doaStd=std(aa);
            if numel(tt)>=2
                d=diff(tt); ds=sort(d); q=max(1,ceil(0.65*numel(ds)));
                st(k).PRI=median(ds(1:q));
                mult=max(1,round(d/max(st(k).PRI,eps)));
                res=d./mult;
                st(k).priStd=std(res);
            else
                st(k).PRI=0; st(k).priStd=0;
            end
            st(k).priAgile=st(k).PRI>0 && st(k).priStd/st(k).PRI>0.018;
            st(k).freqAgile=st(k).freqLevels>=3 || st(k).freqSpan>0.32e6;
            if st(k).priAgile && st(k).freqAgile, st(k).agility='PRI + frequency';
            elseif st(k).priAgile, st(k).agility='PRI agile';
            elseif st(k).freqAgile, st(k).agility='Frequency agile';
            else, st(k).agility='Fixed'; end
        end
    end
