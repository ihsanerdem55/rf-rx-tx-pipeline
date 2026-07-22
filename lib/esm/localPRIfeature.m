    function lpri = localPRIfeature(toa,f,doa,pw,E)
        m=numel(toa); lpri=zeros(m,1);
        pmin=0.35*min([E.PRI]); pmax=2.5*max([E.PRI]);
        for i=1:m
            best=inf; val=median([E.PRI]);
            for j=i+1:min(m,i+12)
                dt=toa(j)-toa(i); if dt<pmin || dt>pmax, continue; end
                c=abs(doa(j)-doa(i))/12 + abs(log(pw(j)/pw(i)))/0.45;
                c=c+0.35*min(abs(f(j)-f(i))/0.8e6,2);
                if c<best, best=c; val=dt; end
            end
            lpri(i)=val;
        end
    end
