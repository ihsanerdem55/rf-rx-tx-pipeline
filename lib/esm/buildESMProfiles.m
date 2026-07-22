    function E = buildESMProfiles(P)
        % Scenario presets deliberately span frequency, PRI, PW, DOA, range and
        % agility. Emitter 1 inherits the normal GUI controls so changes to the
        % original source remain visible in the ESM scene.
        n=P.esmN;
        T=struct('id',0,'name','','classTruth','','fc',0,'rf',0,'PRI',0, ...
            'PW',0,'az',0,'range',0,'amp',1,'modtype','CW','bchirp',0, ...
            'primode','Constant','pridev',0,'fhmode','Fixed','hopstep',0, ...
            'vrad',0);
        E=repmat(T,1,n);

        bw=max(0.5e6,P.Bhi-P.Blo);
        if P.iq
            room=max(0.3e6,min(0.4*P.Fs,0.42*bw));
        else
            room=max(0.3e6,min(P.Fctr-P.Blo,P.Bhi-P.Fctr)-0.55*max(P.bchirp,0));
        end
        df=max(0.15e6,min(1.0e6,0.42*room));
        basePRI=max(P.PRI,25e-6); basePW=min(P.PW,0.45*basePRI);
        azSet=[P.Az, -32, 38, 6];
        rSet=[P.Rng, 24e3, 8e3, 31e3];

        % Primary emitter.
        E(1)=T; E(1).id=1; E(1).name='Primary'; E(1).classTruth='Track radar';
        E(1).fc=P.Fctr; E(1).rf=P.Frf; E(1).PRI=basePRI; E(1).PW=basePW;
        E(1).az=azSet(1); E(1).range=rSet(1); E(1).amp=1;
        E(1).modtype=P.modtype; E(1).bchirp=P.bchirp;
        E(1).primode=P.primode; E(1).pridev=P.pridev;
        E(1).fhmode=P.fhmode; E(1).hopstep=P.hopstep; E(1).vrad=P.vrad;

        for k=2:n
            E(k)=T; E(k).id=k; E(k).name=sprintf('Emitter %d',k);
            E(k).az=azSet(k); E(k).range=rSet(k);
        end

        switch P.esmScenario
            case 'Co-channel MUSIC'
                if n>=2
                    E(2).classTruth='Track radar'; E(2).fc=P.Fctr;
                    E(2).rf=P.Frf; E(2).PRI=0.78*basePRI; E(2).PW=0.55*basePW;
                    E(2).az=max(-75,min(75,P.Az+14)); E(2).amp=0.85;
                    E(2).modtype='CW'; E(2).primode='Constant';
                end
                if n>=3
                    E(3).classTruth='Missile seeker'; E(3).fc=P.Fctr-df;
                    E(3).rf=P.Frf-df; E(3).PRI=max(32e-6,0.48*basePRI);
                    E(3).PW=max(1.2e-6,0.25*basePW); E(3).amp=0.95;
                    E(3).modtype='CW'; E(3).primode='Constant';
                end
                if n>=4
                    E(4).classTruth='Search radar'; E(4).fc=P.Fctr+df;
                    E(4).rf=P.Frf+df; E(4).PRI=1.7*basePRI;
                    E(4).PW=max(4e-6,0.75*basePW); E(4).amp=0.72;
                    E(4).modtype='LFM chirp'; E(4).bchirp=min(0.8e6,max(0.2e6,P.bchirp));
                end
            case 'Dense agile'
                if n>=2
                    E(2).classTruth='Agile multifunction'; E(2).fc=P.Fctr+0.55*df;
                    E(2).rf=P.Frf+0.55*df; E(2).PRI=0.85*basePRI;
                    E(2).PW=max(2e-6,0.45*basePW); E(2).amp=0.9;
                    E(2).primode='Random jitter'; E(2).pridev=0.09;
                    E(2).fhmode='Random hop'; E(2).hopstep=0.35*df;
                end
                if n>=3
                    E(3).classTruth='Missile seeker'; E(3).fc=P.Fctr-0.7*df;
                    E(3).rf=P.Frf-0.7*df; E(3).PRI=max(30e-6,0.5*basePRI);
                    E(3).PW=max(1.0e-6,0.22*basePW); E(3).amp=1.0;
                    E(3).primode='2-level stagger'; E(3).pridev=0.07;
                    E(3).fhmode='3-tone hop'; E(3).hopstep=0.28*df;
                end
                if n>=4
                    E(4).classTruth='Search radar'; E(4).fc=P.Fctr+1.2*df;
                    E(4).rf=P.Frf+1.2*df; E(4).PRI=1.55*basePRI;
                    E(4).PW=max(5e-6,0.75*basePW); E(4).amp=0.68;
                    E(4).modtype='LFM chirp'; E(4).bchirp=min(1.0e6,max(0.3e6,P.bchirp));
                    E(4).primode='Random jitter'; E(4).pridev=0.04;
                end
            otherwise % Mixed air picture
                if n>=2
                    E(2).classTruth='Search radar'; E(2).fc=P.Fctr+df;
                    E(2).rf=P.Frf+df; E(2).PRI=1.65*basePRI;
                    E(2).PW=max(5e-6,0.70*basePW); E(2).amp=0.72;
                    E(2).modtype='LFM chirp'; E(2).bchirp=min(0.9e6,max(0.25e6,P.bchirp));
                    E(2).primode='Constant';
                end
                if n>=3
                    E(3).classTruth='Missile seeker'; E(3).fc=P.Fctr-df;
                    E(3).rf=P.Frf-df; E(3).PRI=max(30e-6,0.52*basePRI);
                    E(3).PW=max(1.1e-6,0.24*basePW); E(3).amp=1.0;
                    E(3).primode='Constant';
                end
                if n>=4
                    E(4).classTruth='Agile multifunction'; E(4).fc=P.Fctr+0.18*df;
                    E(4).rf=P.Frf+0.18*df; E(4).PRI=0.92*basePRI;
                    E(4).PW=max(2.2e-6,0.42*basePW); E(4).amp=0.82;
                    E(4).primode='2-level stagger'; E(4).pridev=0.08;
                    E(4).fhmode='3-tone hop'; E(4).hopstep=0.28*df;
                end
        end

        % Final physical/sampling guards for every preset.
        for k=1:n
            E(k).PRI=max(20e-6,E(k).PRI);
            E(k).PW=max(0.6e-6,min(E(k).PW,0.45*E(k).PRI));
            E(k).az=max(-80,min(80,E(k).az)); E(k).range=max(1e3,E(k).range);
            if P.iq
                lim=0.42*P.Fs;
                E(k).fc=max(-lim,min(lim,E(k).fc));
            else
                lo=P.Blo+0.55*E(k).bchirp; hi=P.Bhi-0.55*E(k).bchirp;
                if hi<=lo, lo=0.05*P.Fs; hi=0.45*P.Fs; end
                E(k).fc=max(lo,min(hi,E(k).fc));
            end
            E(k).rf=P.Flo+E(k).fc;
        end
    end
