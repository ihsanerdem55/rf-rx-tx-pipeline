    function [P, notes] = validate(P)
        % Guard against degenerate / aliasing settings. Rather than crashing or
        % silently producing nonsense, clamp to sane values and record a
        % human-readable note that gets shown in the results panel.
        notes = {};

        if ~(isfinite(P.Fs) && P.Fs > 0)
            P.Fs = 20e6;
            notes{end+1} = 'Sample rate invalid -> reset to 20 MHz.';
        end
        nyq = P.Fs/2;
        if ~(isfinite(P.PRI) && P.PRI > 0)
            P.PRI = 10e-6;
            notes{end+1} = 'PRI invalid -> reset to 10 us.';
        end
        if ~(isfinite(P.PW) && P.PW > 0)
            P.PW = 1e-6;
            notes{end+1} = 'PW invalid -> reset to 1 us.';
        end
        if P.PW >= P.PRI
            P.PW = 0.9*P.PRI;
            notes{end+1} = 'PW >= PRI -> clamped to 90% of PRI.';
        end
        if P.Np < 2
            P.Np = 2;
            notes{end+1} = 'Pulse count raised to 2 (minimum for Doppler processing).';
        elseif P.Np > 256
            P.Np = 256;
            notes{end+1} = 'Pulse count capped at 256 (memory/runtime guard).';
        end

        % --- search band must be ordered and inside [0, Nyquist] ----------
        if P.Bhi < P.Blo
            [P.Blo, P.Bhi] = deal(P.Bhi, P.Blo);
            notes{end+1} = 'Band low/high were reversed -> swapped.';
        end
        if P.Blo < 0,   P.Blo = 0;   notes{end+1} = 'Band low clamped to 0.'; end
        if P.Bhi > nyq, P.Bhi = nyq; notes{end+1} = 'Band high clamped to Nyquist.'; end
        if P.Bhi <= P.Blo
            P.Blo = 0; P.Bhi = nyq;
            notes{end+1} = 'Empty band -> reset to full [0, Nyquist].';
        end

        % --- realisable IIR/FIR guards only apply to a REAL IF bandpass ----
        % (in baseband the receiver filter is an FFT low-pass, see applyFilter).
        if ~P.iq
            switch P.filt
                case {'Butterworth IIR', 'Chebyshev I IIR'}
                    lo = 0.005*P.Fs; hi = 0.495*P.Fs;
                    if P.Blo < lo
                        P.Blo = lo;
                        notes{end+1} = 'IIR: band low raised off DC (design limit).';
                    end
                    if P.Bhi > hi
                        P.Bhi = hi;
                        notes{end+1} = 'IIR: band high lowered off Nyquist (design limit).';
                    end
                    if P.ord > 8
                        P.ord = 8;
                        notes{end+1} = 'IIR order capped at 8 (b/a form gets ill-conditioned).';
                    end
                    P.ord = max(1, P.ord);
                case {'FIR windowed-sinc', 'FIR Kaiser'}
                    P.ord = min(2047, max(15, P.ord));
                    if mod(P.ord, 2) == 0, P.ord = P.ord + 1; end
            end
        else
            P.ord = max(1, P.ord);
        end
        P.tap = min(0.9, max(0.001, P.tap));

        % --- Welch estimator parameters -----------------------------------
        P.wseg = min(64, max(1, P.wseg));      % 1 segment == plain periodogram
        P.wov  = min(0.95, max(0, P.wov));     % keep some non-overlap
        P.crip = min(10, max(0.05, P.crip));   % Chebyshev ripple 0.05..10 dB
        P.kbeta = min(20, max(0, P.kbeta));    % Kaiser beta 0..20

        % --- LFM chirp: keep the sweep inside the sampled band ------------
        if isChirp(P)
            if P.bchirp < 0
                P.bchirp = abs(P.bchirp);
                notes{end+1} = 'Negative chirp BW interpreted as a positive bandwidth.';
            end
            if P.iq
                maxB = 0.9 * P.Fs;             % baseband: |f| <= B/2 < Fs/2
            else
                maxB = 2*min(P.Fif, nyq - P.Fif);  % IF: fc +/- B/2 in (0,nyq)
            end
            maxB = max(0, maxB);
            if P.bchirp > maxB
                P.bchirp = maxB;
                notes{end+1} = sprintf('Chirp BW clamped to %.2f MHz (stays sampled).', maxB/1e6);
            end
            if P.bchirp > (P.Bhi - P.Blo)
                notes{end+1} = 'Chirp BW exceeds the passband -> the filter will clip the sweep.';
            end
            if P.bchirp > 0 && P.bchirp*P.PW < 2
                notes{end+1} = 'LFM time-bandwidth product < 2: pulse-compression gain will be small.';
            end
        end

        chirpHalf = 0;
        if isCompressible(P), chirpHalf = effBW(P)/2; end

        % --- moving emitter, multipath and agility -----------------------
        if ~isfinite(P.vrad), P.vrad = 0; end
        if abs(P.vrad) > 3000
            P.vrad = sign(P.vrad)*3000;
            notes{end+1} = 'Emitter radial velocity capped at +/-3000 m/s.';
        end
        P.fdRx = P.vrad * P.Frf / P.c;         % one-way emitter Doppler (Hz)

        P.mploss = min(80, max(0, P.mploss));
        P.mpamp = 10^(-P.mploss/20);
        if ~(P.mpdelay >= 0), P.mpdelay = 0; end
        if P.mpdelay > 2*P.PRI
            P.mpdelay = 2*P.PRI;
            notes{end+1} = 'Multipath excess delay capped at 2 PRI.';
        end
        P.mpaz = min(90, max(-90, P.mpaz));
        if P.multipath && P.mpdelay == 0
            notes{end+1} = 'Multipath delay is zero: reflection only changes amplitude/DF phase, not range.';
        end
        if P.multipath && isCompressible(P) && effBW(P) > 0 && ...
                P.mpdelay < 1/effBW(P)
            notes{end+1} = 'Multipath delay is below approximately 1/B: direct and reflected range peaks are unresolved.';
        end

        % --- antenna scan modulation -------------------------------------
        if ~isfield(P,'scanmode'), P.scanmode = 'Off'; end
        P.scanrate = max(0, P.scanrate);
        P.scanbw   = min(180, max(0.5, P.scanbw));
        if ~strcmp(P.scanmode,'Off')
            prfTrue = 1/P.PRI;
            if P.scanrate > 0.5*prfTrue
                notes{end+1} = sprintf('Scan rate above PRF/2 (%.0f Hz): the scan-rate estimate will alias.', 0.5*prfTrue);
            end
            if P.scanrate*P.Np*P.PRI < 1
                notes{end+1} = 'Under one scan cycle in the CPI: raise N pulses (or scan rate) to resolve the scan.';
            end
        end

        P.pridev = min(0.45, max(0, P.pridev));
        if strcmp(P.primode,'Constant'), P.pridev = 0; end
        P.hopstep = abs(P.hopstep);
        switch P.fhmode
            case '3-tone hop', hopLevels = 1;
            case 'Random hop', hopLevels = 2;
            otherwise, hopLevels = 0; P.hopstep = 0;
        end
        maxHop = hopLevels*P.hopstep;
        if P.PRI > 0
            minPRI = P.PRI*(1-P.pridev);
            if P.PW >= minPRI
                P.PW = 0.8*minPRI;
                notes{end+1} = 'PW shortened to 80% of the minimum agile PRI.';
            end
        end

        % Keep every hop, receive-Doppler shift and chirp edge inside the ADC.
        if P.iq
            room = max(0, 0.95*nyq - chirpHalf - abs(P.fdRx));
        else
            room = max(0, min(P.Fif-0.005*P.Fs, 0.495*P.Fs-P.Fif) ...
                - chirpHalf - abs(P.fdRx));
        end
        if maxHop > room && hopLevels > 0
            P.hopstep = room/hopLevels;
            maxHop = room;
            notes{end+1} = sprintf('Hop step clamped to %.3f MHz to remain sampled.', P.hopstep/1e6);
        end
        if ~P.iq && (P.Fif-maxHop+P.fdRx < P.Blo || P.Fif+maxHop+P.fdRx > P.Bhi)
            notes{end+1} = 'Some hopped pulses fall outside the receiver search band and may be attenuated/missed.';
        end
        if P.iq && maxHop+abs(P.fdRx)+chirpHalf > 0.5*(P.Bhi-P.Blo)
            notes{end+1} = 'Some IQ hop energy lies outside the configured baseband search width.';
        end
        lambdaMotion = P.c/max(P.Frf,eps);
        vUnamb = lambdaMotion/(2*P.PRI);
        if abs(P.vrad) > vUnamb
            notes{end+1} = sprintf('Receive Doppler aliases at the nominal PRF (unambiguous velocity +/-%.1f m/s).', vUnamb);
        end
        if ~strcmp(P.fhmode,'Fixed')
            notes{end+1} = 'With frequency hopping, the global carrier peak represents one channel; use the per-pulse tracker in tab 11.';
        end

        % --- multi-emitter ESM --------------------------------------------
        P.esmN = min(4,max(2,P.esmN));
        % Dedicated ESM array geometry: a half-wavelength uniform linear
        % array at the nominal RF. MUSIC needs ~lambda/2 spacing; reusing the
        % interferometer's sparse fine baseline (>lambda) causes grating lobes.
        P.esmPos = [0, 0.5*P.c/max(P.Frf,eps), P.c/max(P.Frf,eps)];
        P.musicStep = min(2,max(0.1,P.musicStep));
        P.pdwMiss = min(0.45,max(0,P.pdwMiss));
        if P.esmOn && P.esmN > 2
            notes{end+1} = 'Three sensors support at most two simultaneous co-channel MUSIC sources; additional emitters are separated by frequency/PDW tracks.';
        end

        % --- DRFM deception / EW duel -------------------------------------
        P.ntgt    = min(12, max(1, P.ntgt));
        if ~(P.drdelay >= 0), P.drdelay = 0; end
        if abs(P.dop) >= nyq, P.dop = 0; notes{end+1} = 'Doppler magnitude too large -> set to 0.'; end
        if ~isfinite(P.vgstep), P.vgstep = 0; end
        maxVGstep=0.45*P.Fs/max(1,P.Np-1);
        P.vgstep = min(maxVGstep,max(-maxVGstep,P.vgstep));
        P.coverN = min(20,max(1,P.coverN));
        P.coverLoss = min(60,max(0,P.coverLoss));
        P.coverAmp = 10^(-P.coverLoss/20);
        P.jsdb = min(60,max(-30,P.jsdb));
        P.spotbw = min(0.9*P.Fs,max(P.Fs/4096,abs(P.spotbw)));
        P.vpriDev = min(0.45,max(0,P.vpriDev));
        P.vhopstep = min(max(0,0.45*P.Fs-chirpHalf),max(0,abs(P.vhopstep)));
        P.clutterDB = min(50,max(-30,P.clutterDB));
        P.clutterV = min(100,max(0,abs(P.clutterV)));
        P.pfa = min(0.1,max(1e-8,P.pfa));
        if P.eccmPRI && P.vpriDev == 0
            notes{end+1} = 'Victim PRI ECCM selected with zero stagger depth.';
        end
        if P.eccmFH && P.vhopstep == 0
            notes{end+1} = 'Victim frequency-agility ECCM selected with zero hop step.';
        end
        fdEnd=P.dop+(P.Np-1)*P.vgstep;
        if max(abs([P.dop fdEnd]))>=nyq
            notes{end+1}='VGPO endpoint exceeds the sampled fast-time Nyquist interval and will alias.';
        end
        if P.cover && ~strcmp(P.txmode,'DRFM (capture+shift)')
            notes{end+1} = 'Cover pulses are applied only in DRFM mode.';
        end
        if ~strcmp(P.jammode,'Off') && ~strcmp(P.txmode,'DRFM (capture+shift)')
            notes{end+1} = 'Noise jamming is applied only in DRFM mode.';
        end
        if ~strcmp(P.mtiMode,'Off') && P.Np < 4
            notes{end+1} = 'MTI has very little slow-time support with fewer than four pulses.';
        end

        % --- IF -> RF upconversion sanity ---------------------------------
        if ~(P.Ftrf > 0)
            P.Ftrf = 2e9; notes{end+1} = 'Tx RF invalid -> set to 2 GHz.';
        end
        P.FloTx = P.Ftrf - P.Ft;               % recompute LO after any clamp
        if P.Ftrf < P.Ft
            notes{end+1} = 'Tx RF below Tx IF (upconversion LO is negative).';
        end

        % --- RF front-end / domain sanity ---------------------------------
        chirpHalf = 0; if isCompressible(P), chirpHalf = effBW(P)/2; end
        if P.iq
            % Baseband: content spans +/-(chirpHalf); need that < Nyquist.
            if chirpHalf >= nyq
                notes{end+1} = 'Baseband: signal wider than Fs -> aliasing (raise Fs).';
            end
        else
            % IF: content spans Fif +/- chirpHalf; must stay in (0, Nyquist).
            if P.Fif <= 0
                notes{end+1} = 'IF <= 0: pick Baseband IQ for a zero-IF receiver.';
            end
            if P.Fif + chirpHalf >= nyq
                notes{end+1} = sprintf(['IF+chirp (%.2f MHz) >= Nyquist %.2f MHz: ' ...
                    'aliasing. Raise Fs above %.1f MHz.'], ...
                    (P.Fif+chirpHalf)/1e6, nyq/1e6, 2*(P.Fif+chirpHalf)/1e6);
            end
            if P.Fif < P.Blo || P.Fif > P.Bhi
                notes{end+1} = 'IF carrier is OUTSIDE the search band (Rx will miss it).';
            end
        end
        if P.Frf <= P.Fctr
            notes{end+1} = 'RF carrier must exceed the IF (downconversion LO <= 0).';
        end
        if P.Ft >= nyq
            notes{end+1} = sprintf('Tx freq %.3g MHz >= Nyquist %.3g MHz (aliasing).', ...
                P.Ft/1e6, nyq/1e6);
        end
        if ~(P.PW > 0),  P.PW  = 1e-6;  notes{end+1} = 'PW <= 0 -> set to 1 us.'; end
        if ~(P.PRI > 0), P.PRI = 10e-6; notes{end+1} = 'PRI <= 0 -> set to 10 us.'; end
        if P.PW >= P.PRI
            P.PW = 0.9 * P.PRI;
            notes{end+1} = 'PW >= PRI -> clamped to 90% of PRI.';
        end
        % Prevent accidental multi-million-sample DRFM buffers while still
        % allowing a large pull-off spanning several unambiguous range windows.
        if P.drdelay > 4*P.PRI
            P.drdelay = 4*P.PRI;
            notes{end+1} = 'RGPO delay capped at 4 PRI (memory/runtime guard).';
        end

        % --- source geometry / array sanity -------------------------------
        if ~(P.Rng >= 0), P.Rng = 0; notes{end+1} = 'Range < 0 -> set to 0.'; end
        if P.Az >  90, P.Az =  90; notes{end+1} = 'Azimuth clamped to +90 deg (array limit).'; end
        if P.Az < -90, P.Az = -90; notes{end+1} = 'Azimuth clamped to -90 deg (array limit).'; end
        if ~(P.base > 0)
            P.base = 0.1; notes{end+1} = 'Fine baseline <= 0 -> set to 0.1 m.';
        end
        if ~(P.bshort > 0), P.bshort = 0.01; end
        P.phoff  = min(pi, max(0, P.phoff));
        P.gainmm = min(0.9, max(0, P.gainmm));
        % DF ambiguity uses the true RF wavelength (physical, not the IF).
        lambda = P.c / max(P.Frf, eps);
        dual = strcmp(P.dfmode, 'Dual-baseline');
        if dual
            % The SHORT baseline must stay <= lambda/2 or it cannot resolve the
            % fine baseline's wrap (its own coarse angle would alias).
            if P.bshort > lambda/2
                notes{end+1} = sprintf(['Short baseline %.3g m > lambda/2 = %.3g m: ' ...
                    'ambiguity resolution may fail (shorten it).'], P.bshort, lambda/2);
            end
            if P.base <= P.bshort
                notes{end+1} = 'Fine baseline <= short baseline: no resolution benefit (make fine larger).';
            end
        else
            % Single baseline: honest warning that the estimate wraps beyond
            % +/-asin(lambda/2/base) when the baseline exceeds lambda/2.
            if P.base > lambda/2
                notes{end+1} = sprintf(['Single-baseline & base %.3g m > lambda/2 = %.3g m: ' ...
                    'azimuth AMBIGUOUS beyond +/-%.0f deg (use Dual-baseline).'], ...
                    P.base, lambda/2, asind(min(1, lambda/(2*P.base))));
            end
        end
        if P.Rng / P.c >= P.PRI
            notes{end+1} = 'Range delay >= PRI: unambiguous here only by the synchronised-time assumption.';
        end
    end
