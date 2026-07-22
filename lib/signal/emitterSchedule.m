    function [emitTimes, priTrue, hopTrue] = emitterSchedule(P, seed)
        % Deterministic schedule for a given seed. PRI agility changes pulse
        % start times; frequency agility changes the carrier independently on
        % each pulse. The patterns are zero-mean so the entered PRI/RF remain
        % the nominal centre values.
        rs = RandStream('mt19937ar','Seed',seed+4242);
        priTrue = P.PRI*ones(1,max(0,P.Np-1));
        switch P.primode
            case '2-level stagger'
                pat = repmat([-1 1],1,ceil(numel(priTrue)/2));
                priTrue = P.PRI*(1+P.pridev*pat(1:numel(priTrue)));
            case 'Random jitter'
                priTrue = P.PRI*(1+P.pridev*(2*rand(rs,1,numel(priTrue))-1));
        end
        emitTimes = [0 cumsum(priTrue)];

        switch P.fhmode
            case '3-tone hop'
                lev = mod(0:P.Np-1,3)-1;      % -1, 0, +1, ...
            case 'Random hop'
                lev = floor(5*rand(rs,1,P.Np))-2; % integer channels -2..+2
            otherwise
                lev = zeros(1,P.Np);
        end
        hopTrue = lev*P.hopstep;
    end
