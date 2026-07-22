    function [t, gate, clean, rx1, rx2, rx3, M] = makeSignal(P, seed)
        % Build a pulsed moving emitter with optional PRI agility, frequency
        % hopping and a delayed ground-reflection path. The RF itself is still
        % represented by a sampled IF/baseband-equivalent waveform, while the
        % true RF frequency determines Doppler and array phase.
        Fs = P.Fs;
        [emitTimes, priTrue, hopTrue] = emitterSchedule(P, seed);

        % Positive radial velocity means closing, so the source range decreases
        % slightly across the coherent processing interval. Its one-way Doppler
        % is fd=v/lambda, not the two-way 2v/lambda used by monostatic radar.
        ranges = max(0, P.Rng - P.vrad*emitTimes);
        tau = ranges/P.c;
        arrival = emitTimes + tau;
        mpArrival = arrival + P.mpdelay;
        tailPRI = P.PRI;
        if ~isempty(priTrue), tailPRI = max(priTrue); end
        lastPath = arrival(end);
        if P.multipath, lastPath = mpArrival(end); end
        Ntot = round((lastPath + P.PW + 0.5*tailPRI)*Fs) + 2;
        if mod(Ntot,2) == 1, Ntot = Ntot + 1; end
        t = (0:Ntot-1)/Fs;

        % Fixed hardware errors are generated from a hardware-only seed so an
        % array calibration remains valid when the user requests new noise.
        hw = RandStream('mt19937ar', 'Seed', 90210);
        phE = zeros(1,3); gnE = ones(1,3);
        if P.phoff  > 0, phE(2:3) = P.phoff*randn(hw,1,2); end
        if P.gainmm > 0, gnE(2:3) = 1 + P.gainmm*randn(hw,1,2); end

        gate = zeros(1,Ntot);                   % direct-path truth gate
        c1 = complex(zeros(1,Ntot));
        c2 = complex(zeros(1,Ntot));
        c3 = complex(zeros(1,Ntot));
        lfm = strcmp(P.modtype,'LFM chirp');
        K = 0; if lfm, K = P.bchirp/P.PW; end

        fcProcTrue = P.Fctr + hopTrue + P.fdRx;
        rfTrue = P.Frf + hopTrue + P.fdRx;
        for k = 1:P.Np
            fproc = fcProcTrue(k);
            lambda = P.c/max(rfTrue(k),eps);

            % ---------------- direct path ----------------
            t0 = arrival(k);
            idx = find(t >= t0 & t < t0+P.PW);
            if ~isempty(idx)
                gate(idx) = 1;
                tl = t(idx)-(t0+P.PW/2);
                intra = zeros(size(tl));
                if lfm, intra = pi*K*tl.^2; end
                b = exp(1i*(2*pi*fproc*t(idx)+intra));
                dF = 2*pi*P.base*sind(P.Az)/lambda;
                dS = 2*pi*P.bshort*sind(P.Az)/lambda;
                c1(idx) = c1(idx)+b;
                c2(idx) = c2(idx)+gnE(2)*b.*exp(1i*(dF+phE(2)));
                c3(idx) = c3(idx)+gnE(3)*b.*exp(1i*(dS+phE(3)));
            end

            % ---------------- reflected path ----------------
            if P.multipath && P.mpamp > 0
                tm = mpArrival(k);
                im = find(t >= tm & t < tm+P.PW);
                if ~isempty(im)
                    tlm = t(im)-(tm+P.PW/2);
                    intraM = zeros(size(tlm));
                    if lfm, intraM = pi*K*tlm.^2; end
                    % Ground reflection: excess propagation phase plus a pi
                    % phase reversal. Using the RF phase makes small path-delay
                    % changes produce realistic constructive/destructive fading.
                    reflPhase = -2*pi*rfTrue(k)*P.mpdelay + pi;
                    bm = P.mpamp*exp(1i*(2*pi*fproc*t(im)+intraM+reflPhase));
                    dFm = 2*pi*P.base*sind(P.mpaz)/lambda;
                    dSm = 2*pi*P.bshort*sind(P.mpaz)/lambda;
                    c1(im) = c1(im)+bm;
                    c2(im) = c2(im)+gnE(2)*bm.*exp(1i*(dFm+phE(2)));
                    c3(im) = c3(im)+gnE(3)*bm.*exp(1i*(dSm+phE(3)));
                end
            end
        end

        if P.iq
            clean = c1; clean2 = c2; clean3 = c3;
        else
            clean = real(c1); clean2 = real(c2); clean3 = real(c3);
        end

        on = clean(gate > 0);
        if isempty(on), Ps = 0.5; else, Ps = mean(abs(on).^2); end
        Pn = Ps/10^(P.SNR/10);
        s1 = RandStream('mt19937ar','Seed',seed);
        s2 = RandStream('mt19937ar','Seed',seed+777);
        s3 = RandStream('mt19937ar','Seed',seed+1554);
        if P.iq
            rx1 = clean  + sqrt(Pn/2)*(randn(s1,1,Ntot)+1i*randn(s1,1,Ntot));
            rx2 = clean2 + sqrt(Pn/2)*(randn(s2,1,Ntot)+1i*randn(s2,1,Ntot));
            rx3 = clean3 + sqrt(Pn/2)*(randn(s3,1,Ntot)+1i*randn(s3,1,Ntot));
        else
            rx1 = clean  + sqrt(Pn)*randn(s1,1,Ntot);
            rx2 = clean2 + sqrt(Pn)*randn(s2,1,Ntot);
            rx3 = clean3 + sqrt(Pn)*randn(s3,1,Ntot);
        end

        M.emitTimes = emitTimes;
        M.priTrue = priTrue;
        M.hopTrue = hopTrue;
        M.fcProcTrue = fcProcTrue;
        M.rfTrue = rfTrue;
        M.ranges = ranges;
        M.arrivalTimes = arrival;
        M.mpArrivalTimes = mpArrival;
    end
