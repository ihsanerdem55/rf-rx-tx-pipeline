    function [emit,pris,hops] = esmEmitterSchedule(E,Tobs,rs)
        emit=0; pris=[]; hops=[]; k=1;
        while emit(end)<Tobs && k<=96
            switch E.fhmode
                case '3-tone hop', lev=mod(k-1,3)-1;
                case 'Random hop', lev=floor(5*rand(rs))-2;
                otherwise, lev=0;
            end
            hops(k)=lev*E.hopstep; %#ok<AGROW>
            switch E.primode
                case '2-level stagger', d=E.PRI*(1+E.pridev*(-1)^k);
                case 'Random jitter', d=E.PRI*(1+E.pridev*(2*rand(rs)-1));
                otherwise, d=E.PRI;
            end
            pris(k)=d; %#ok<AGROW>
            emit(k+1)=emit(k)+d; %#ok<AGROW>
            k=k+1;
        end
        emit=emit(emit<Tobs); hops=hops(1:numel(emit));
        pris=pris(1:max(0,numel(emit)-1));
    end
