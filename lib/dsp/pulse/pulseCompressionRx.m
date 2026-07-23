    function C = pulseCompressionRx(z, fcCtr, P, riseHint, fcPulse)
        % Match-filter the received LFM pulse train with the ideal reference
        % chirp. Everything is implemented with base MATLAB: complex mixing,
        % convolution, local peak searches and a three-point peak refinement.
        %
        % The correlation output time is registered to the START of the pulse:
        % if a reference of L samples starts at input sample n0, its matched-
        % filter peak occurs at output sample n0+L-1. Therefore
        %       t_pc = (output_index - L)/Fs
        % maps that peak back to the original pulse-onset time.
        C = struct('active',false,'t',[],'mag',[],'magDB',[], ...
            'peakIdx',[],'peakTimes',[],'firstPeakIdx',NaN, ...
            'tauEst',NaN,'width3dB',NaN,'TB',NaN,'gainTheory',NaN, ...
            'toaRes',NaN,'radarRes',NaN,'snrIn',NaN,'snrOut',NaN, ...
            'gainMeasured',NaN,'Ncoh',0,'psl',NaN,'pcWin','None');

        if nargin < 4, riseHint = []; end
        if nargin < 5, fcPulse = []; end
        Beff = effBW(P);                      % sweep BW (chirp) or chip rate (code)
        if ~isCompressible(P) || ~(Beff > 0) || isempty(z)
            return
        end

        Fs = P.Fs;
        L = max(2,round(P.PW*Fs));
        tp = (0:L-1)/Fs;
        tl = tp - P.PW/2;
        ref = exp(1i*chirpPhase(tl, P));      % unit-amplitude chirp or phase code
        pcwin = 'None'; if isfield(P,'pcwin'), pcwin = P.pcwin; end
        ref = ref .* pcWindow(L, pcwin);      % receive taper for sidelobe control
        C.pcWin = pcwin;
        Eref = sum(abs(ref).^2);

        % Remove carrier and hop frequency. When pulse-wise estimates are
        % available, overwrite each pulse neighbourhood with its own dehopping
        % oscillator; otherwise use the global occupied-band centre.
        if ~isfinite(fcCtr), fcCtr = P.Fctr; end
        nabs = 0:numel(z)-1;
        bb = z(:).' .* exp(-1i*2*pi*fcCtr*nabs/Fs);
        if ~isempty(riseHint) && numel(fcPulse) >= numel(riseHint)
            ext = round((P.PW+P.mpdelay+0.15*P.PRI)*Fs);
            for q = 1:numel(riseHint)
                ia = max(1,riseHint(q)-round(0.1*P.PW*Fs));
                ib = min(numel(z),riseHint(q)+ext);
                nn = ia-1:ib-1;
                bb(ia:ib) = z(ia:ib).*exp(-1i*2*pi*fcPulse(q)*nn/Fs);
            end
        end
        mf = conv(bb,conj(fliplr(ref)),'full');
        mag = abs(mf)/(Eref+eps);
        tpc = ((1:numel(mf))-L)/Fs;

        % With detected pulse hints, search independently around every actual
        % (possibly staggered) arrival. Without hints, retain the coherent
        % constant-PRI scorer that is valuable at very low input SNR.
        if ~isempty(riseHint)
            targets = riseHint+L-1;
            ncoh = numel(targets);
            seed = targets(1);
        else
            step = max(1,round(P.PRI*Fs));
            ncoh = max(1,min(P.Np,floor((numel(mag)-1)/step)+1));
            maxStart = numel(mag)-(ncoh-1)*step;
            score = zeros(1,maxStart);
            for q = 0:ncoh-1
                score = score+mag((1:maxStart)+q*step).^2;
            end
            valid = find(tpc(1:maxStart)>=0);
            if isempty(valid), return; end
            [~,jj] = max(score(valid));
            seed = valid(jj);
            targets = seed+(0:ncoh-1)*step;
        end

        % Refine each expected PRI-spaced peak in a narrow local window, then
        % use three-point parabolic interpolation for fractional-sample TOA.
        localHalf = max(3,min(round(0.2*P.PW*Fs), ...
            max(3,round(3*Fs/Beff))));
        targets = targets(targets>=1 & targets<=numel(mag));
        nPeak = numel(targets);
        pidx = zeros(1,nPeak); ptime = NaN(1,nPeak);
        for q = 1:nPeak
            target = targets(q);
            ia = max(1,target-localHalf);
            ib = min(numel(mag),target+localHalf);
            [~,jr] = max(mag(ia:ib));
            ip = ia+jr-1; pidx(q) = ip;
            delta = 0;
            if ip > 1 && ip < numel(mag)
                a = mag(ip-1)^2; b = mag(ip)^2; c = mag(ip+1)^2;
                den = a-2*b+c;
                if abs(den) > eps
                    delta = max(-0.5,min(0.5,0.5*(a-c)/den));
                end
            end
            ptime(q) = (ip+delta-L)/Fs;
        end
        if isempty(pidx) || ~isfinite(ptime(1)), return; end

        ip = pidx(1);
        peak = mag(ip);
        thr3 = peak/sqrt(2);                   % amplitude -3 dB threshold
        il = ip; while il > 1 && mag(il-1) >= thr3, il = il-1; end
        ir = ip; while ir < numel(mag) && mag(ir+1) >= thr3, ir = ir+1; end

        % Peak sidelobe level (PSL): the highest response BEYOND the mainlobe
        % nulls (not the -3 dB edges -- the skirt between the -3 dB points and
        % the first null is part of the mainlobe and would swamp the true
        % sidelobes). Walk out from the peak until the magnitude stops falling
        % to locate the first null on each side, then search a local window so
        % the PRI-spaced neighbouring compressed peaks are not counted.
        nl = ip; while nl > 1        && mag(nl-1) <= mag(nl), nl = nl-1; end
        nr = ip; while nr < numel(mag) && mag(nr+1) <= mag(nr), nr = nr+1; end
        psHalf = max(nr-nl+1, round(min(0.4*P.PRI, 25/max(Beff,eps))*Fs));
        ia2 = max(1,ip-psHalf); ib2 = min(numel(mag),ip+psHalf);
        reg = ia2:ib2; reg = reg(reg < nl | reg > nr);
        if isempty(reg)
            C.psl = NaN;
        else
            C.psl = 20*log10(max(mag(reg))/max(peak,eps) + eps);
        end

        % Estimate input/output SNR using protected pulse/peak regions. These are
        % display diagnostics; the theoretical coherent gain remains 10log10(BT).
        envBB = abs(bb);
        sigMaskIn = false(size(envBB));
        for q = 1:numel(ptime)
            ia = max(1,round(ptime(q)*Fs)+1);
            ib = min(numel(envBB),ia+L-1);
            if ia <= ib, sigMaskIn(ia:ib) = true; end
        end
        noiseIn = envBB(~sigMaskIn);
        if isempty(noiseIn), noiseIn = envBB; end
        noiseAmpIn = median(noiseIn)/sqrt(log(2)+eps); % Rayleigh median -> scale
        if any(sigMaskIn), sigAmpIn = median(envBB(sigMaskIn)); else, sigAmpIn = max(envBB); end
        C.snrIn = 20*log10(max(sigAmpIn,eps)/max(noiseAmpIn,eps));

        sigMaskOut = false(size(mag));
        protect = max(2,round(2*Fs/Beff));
        for q = 1:numel(pidx)
            ia = max(1,pidx(q)-protect); ib = min(numel(mag),pidx(q)+protect);
            sigMaskOut(ia:ib) = true;
        end
        noiseOut = mag(~sigMaskOut);
        if isempty(noiseOut), noiseOut = mag; end
        noiseAmpOut = median(noiseOut)/sqrt(log(2)+eps);
        C.snrOut = 20*log10(max(peak,eps)/max(noiseAmpOut,eps));
        C.gainMeasured = C.snrOut-C.snrIn;
        C.Ncoh = ncoh;

        C.active = true;
        C.t = tpc;
        C.mag = mag;
        C.magDB = 20*log10(mag/(max(mag)+eps)+eps);
        C.magDB = max(C.magDB,-70);
        C.peakIdx = pidx;
        C.peakTimes = ptime;
        C.firstPeakIdx = ip;
        C.tauEst = ptime(1);
        C.width3dB = (ir-il+1)/Fs;
        C.TB = P.PW*Beff;
        C.gainTheory = 10*log10(max(C.TB,1));
        C.toaRes = P.c/Beff;               % one-way R=c*tau resolution scale
        C.radarRes = P.c/(2*Beff);         % conventional monostatic radar
    end
