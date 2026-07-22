    function [az, azc, kw, wrapped] = resolveAz(dphiF, dphiS, lambda, P)
        % Turn measured fine/short inter-channel phases into an azimuth.
        %   Dual-baseline: the SHORT pair (<= lambda/2) gives a coarse but
        %   unambiguous angle that fixes the wrap integer k of the accurate
        %   FINE pair. Single: invert the fine phase directly (ambiguous if
        %   base > lambda/2).
        azc = NaN; kw = 0;
        twoBase = strcmp(P.dfmode, 'Dual-baseline') && P.bshort > 0;
        if twoBase
            sc = dphiS * lambda / (2*pi*P.bshort);
            azc = asind(max(-1, min(1, sc)));
            dphiF_pred = 2*pi*P.base * sc / lambda;
            kw = round((dphiF_pred - dphiF) / (2*pi));
            sin_th = (dphiF + 2*pi*kw) * lambda / (2*pi*P.base);
            wrapped = abs(sin_th) > 1;
        else
            sin_th = dphiF * lambda / (2*pi*P.base);
            wrapped = (P.base > lambda/2) || (abs(sin_th) > 1);
        end
        az = asind(max(-1, min(1, sin_th)));
    end
