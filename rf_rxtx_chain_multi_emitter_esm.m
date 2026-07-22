function rf_rxtx_chain_multi_emitter_esm
% RF_RXTX_CHAIN  Interactive RF receive / retransmit simulation (educational DSP demo).
%
%   A single-file MATLAB app that walks a pulsed RF signal through a complete
%   pipeline and draws graphics at every stage:
%
%     GENERATE : build a noisy pulse train -- either CW or an intra-pulse LFM
%                chirp -- with a user-set TRUE carrier frequency, PRI, pulse
%                width (PW), pulse count and SNR, emitted from a source at a
%                TRUE (range, azimuth) and selectable radial velocity, with
%                optional multipath, PRI stagger/jitter and pulse-to-pulse
%                carrier hopping. Motion creates genuine receive Doppler and
%                slow-time phase progression. The pulse arrives
%                delayed by tau = range/c and is captured on a 3-element
%                receive array, so the auxiliary channels differ from the
%                reference by phase 2*pi*d*sin(azimuth)/lambda.
%                The true carrier / range / azimuth are UNKNOWN to the receiver.
%
%     Rx       : (1) band-limit to a user-defined SEARCH BAND using a
%                    USER-SELECTABLE filter (see FILTERS below),
%                (2) ESTIMATE the carrier as the strongest in-band spectral
%                    peak -- via either a single windowed periodogram (FFT) or
%                    an averaged Welch PSD (both hand-rolled) -- refined with
%                    parabolic (sub-bin) interpolation,
%                (3) MEASURE pulse width (PW) and pulse repetition interval
%                    (PRI) from the signal envelope (analytic-signal detector),
%                (4) for LFM, MATCH-FILTER the complex envelope against the
%                    reference chirp. Pulse compression turns the wide pulse
%                    into a narrow correlation peak, provides the time-bandwidth
%                    processing gain, and sharpens the time-of-arrival estimate,
%                (5) GEOLOCATE the source: RANGE from the first-pulse TOA
%                    (compressed-peak TOA for LFM, envelope edge for CW), and
%                    AZIMUTH from the measured inter-channel phase difference.
%
%     Tx       : RETRANSMIT, either as a DRFM (Digital RF Memory) that stores
%                the captured samples and replays them frequency-translated to a
%                target frequency (preserving envelope, timing and intra-pulse
%                modulation), with optional deception -- Doppler shift, range
%                gate pull-off (delay) and multiple false targets -- or by
%                PARAMETRIC regeneration of a clean pulse train from PW/PRI.
%                The processed IF signal is then UPCONVERTED to an RF carrier
%                (GHz): a transmit mixer (LO = Ftrf - Tx_IF) places it at the
%                RF frequency. As on the Rx side the GHz waveform is never
%                sampled -- the samples are its baseband-equivalent and the RF
%                carrier is carried as a scalar (keeps Fs small and fast).
%
%     VICTIM   : model a coherent pulse-Doppler victim radar. Its known transmit
%                pulse is matched-filtered in fast time (range compression), and
%                an FFT across pulses forms a range-Doppler map. A weak skin
%                return is plotted with the DRFM false targets at their apparent
%                radar ranges R + c*delay/2 and Doppler-derived radial velocity.
%                Doppler/velocity aliasing at the victim PRF is shown honestly.
%
% -------------------------------------------------------------------------
%   FILTERS (all implemented from scratch -- NO Signal Processing Toolbox)
%
%     'FFT brick-wall'      Zero out every FFT bin outside the band. Perfectly
%                           flat, infinitely sharp, zero-phase -- but a hard
%                           spectral edge is a sinc in time, so pulse edges
%                           ring (Gibbs). Non-causal / not realisable.
%
%     'FFT cosine-taper'    Same idea, but the mask edges are raised-cosine
%                           instead of a cliff. Trades a little selectivity
%                           for markedly less time-domain ringing.
%
%     'FFT Gaussian'        FFT mask with a Gaussian band shape (-3 dB at the
%                           band edges). Infinitely smooth transition, so the
%                           time response has virtually no ringing -- the
%                           opposite extreme from the brick-wall.
%
%     'FIR windowed-sinc'   Real, realisable Type-I linear-phase FIR: the ideal
%                           bandpass impulse response truncated to N taps and
%                           Hamming-windowed. Applied with exact group-delay
%                           compensation, so it is effectively zero-phase.
%                           ~ -53 dB stopband with a Hamming window.
%
%     'FIR Kaiser'          Same FIR structure but a Kaiser window (beta = 8)
%                           for a much cleaner stopband (~ -84 dB) at the cost
%                           of a slightly wider transition band.
%
%     'Butterworth IIR'     Classic maximally-flat IIR, designed here by
%                           analog prototype -> lowpass-to-bandpass transform
%                           -> bilinear transform, then applied FORWARD and
%                           BACKWARD (filtfilt-style) for zero phase. Cheap
%                           (few coefficients) but the double pass squares the
%                           magnitude, so the signal sees |H|^2.
%
%     'Chebyshev I IIR'     Same design path with Chebyshev prototype poles
%                           (1 dB equiripple passband). Sharper roll-off than
%                           Butterworth for the same order, at the price of
%                           passband ripple. Also applied zero-phase (|H|^2).
%
%   Note: sinc, butter, cheby1, freqz, filtfilt, impz, kaiser and hilbert are
%   all Signal Processing Toolbox functions, so every one of them is
%   re-implemented locally here. Only base MATLAB (fft, ifft, filter, conv,
%   poly, besseli, movmean) is used.
%
%   Requires R2021a or later (uifigure, scrollable uigridlayout, uitabgroup).
%
%   USAGE:  save as rf_rxtx_chain_ew_duel.m, then at the MATLAB prompt type:
%               >> rf_rxtx_chain_ew_duel
%
%   Units in the GUI: frequencies in MHz, times in microseconds (us).
% -------------------------------------------------------------------------

% ------------------------------------------------------------------ state
% Add the modular library tree to the path so the
% helper functions in lib/ resolve when the app runs.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'lib')));

S = struct();     % computed signals / results (the "model")
h = struct();     % UI handles                (the "view")
fprintf('Launching rf_rxtx_chain_ew_duel: victim ECCM, clutter/MTI, CFAR, VGPO, cover pulses, and noise jamming enabled.\n');

% Array calibration: the measured inter-channel phase (and gain) offsets that
% a boresight reference reveals, to be subtracted before direction finding.
cal = struct('valid', false, 'dphiF', 0, 'dphiS', 0, 'gF', 1, 'gS', 1);

buildGUI();
runPipeline(false);    % populate everything with the defaults on launch

% =====================================================================
%                              GUI
% =====================================================================
    function buildGUI()
        h.fig = uifigure('Name', ...
            'RF Rx/Tx Chain MULTI-EMITTER ESM  |  MUSIC, Deinterleaving, Classification and EW Duel', ...
            'Position', [60 60 1320 820]);

        outer = uigridlayout(h.fig, [1 2]);
        outer.ColumnWidth = {470, '1x'};
        outer.RowHeight   = {'1x'};

        % ---------------- LEFT: controls ----------------
        % Compact TWO-CONTROLS-PER-ROW layout: a 4-column grid (label,field,
        % label,field) roughly halves the vertical height versus one control
        % per row. The grid itself is Scrollable, and every row uses a FIXED
        % height (no '1x'), so when the total content is taller than the
        % viewport a scrollbar appears -- fixing panels that don't fully fit.
        ctrlPanel = uipanel(outer, 'Title', 'Controls', 'Scrollable', 'on');
        cg = uigridlayout(ctrlPanel, [58 4]);
        cg.Scrollable = 'on';               % <-- scrollbar on the control grid
        % Fixed-height rows keep the expanded EW-duel controls predictable;
        % the parent panel remains scrollable.
        cg.RowHeight = repmat({28}, 1, 58);
        cg.RowHeight{57} = 32;               % run buttons
        cg.RowHeight{58} = 1320;             % results / diagnostics
        cg.ColumnWidth = {90, '1x', 90, '1x'};
        cg.RowSpacing = 4;
        cg.Padding = [8 8 8 8];

        r = 1;
        sec4(cg, r, 'SIGNAL GENERATION  (ground truth)'); r = r+1;
        h.Fs  = numL(cg, r, 'Fs (MHz)',      40);      h.Frf = numR(cg, r, 'RF (GHz)',  2);    r = r+1;
        h.Dom = dropFull(cg, r, 'Domain', {'IF (real)', 'Baseband IQ'}, 'IF (real)'); r = r+1;
        h.Fif = numL(cg, r, 'IF (MHz)',       4);                                              r = r+1;
        h.PRI = numL(cg, r, 'PRI (us)',     100);      h.PW  = numR(cg, r, 'PW (us)',   10);   r = r+1;
        h.Np  = numL(cg, r, 'N pulses',       8, '%d'); h.SNR = numR(cg, r, 'SNR (dB)',  15);  r = r+1;
        h.Mod = dropFull(cg, r, 'Mod type', {'CW', 'LFM chirp'}, 'CW'); r = r+1;
        h.Bch = numL(cg, r, 'Chirp BW MHz',   3);                                              r = r+1;

        sec4(cg, r, 'SOURCE GEOMETRY  (ground truth)'); r = r+1;
        h.Rng = numL(cg, r, 'Range (km)',    15);      h.Az  = numR(cg, r, 'Azim (deg)', 20);  r = r+1;

        sec4(cg, r, 'EMITTER MOTION / MULTIPATH'); r = r+1;
        h.Vrad = numL(cg, r, 'Rad vel m/s',   0);      h.MPloss = numR(cg, r, 'MP loss dB', 10); r = r+1;
        h.MPdly = numL(cg, r, 'MP extra us',  3);      h.MPaz = numR(cg, r, 'MP az deg', -15);  r = r+1;
        h.MPmode = dropFull(cg, r, 'Multipath', {'Off', 'On'}, 'Off'); r = r+1;

        sec4(cg, r, 'EMITTER PRI / FREQUENCY AGILITY'); r = r+1;
        h.PRImode = dropFull(cg, r, 'PRI mode', ...
            {'Constant', '2-level stagger', 'Random jitter'}, 'Constant'); r = r+1;
        h.PRIdev = numL(cg, r, 'PRI dev (%)', 10);     h.Hstep = numR(cg, r, 'Hop step MHz', 0.5); r = r+1;
        h.FHmode = dropFull(cg, r, 'Freq agility', ...
            {'Fixed', '3-tone hop', 'Random hop'}, 'Fixed'); r = r+1;

        sec4(cg, r, 'MULTI-EMITTER ESM  (composite scene + PDW processing)'); r = r+1;
        h.ESMmode = dropFull(cg, r, 'ESM scene', {'Off','On'}, 'Off'); r = r+1;
        h.ESMN = numL(cg, r, 'Emitters', 3, '%d');
        h.ESMscenario = dropFullAt(cg, r, 3, 4, 'Scenario', ...
            {'Mixed air picture','Co-channel MUSIC','Dense agile'}, 'Mixed air picture'); r = r+1;
        h.MusicStep = numL(cg, r, 'MUSIC step deg', 0.25);
        h.PDWMiss = numR(cg, r, 'PDW miss (%)', 4); r = r+1;

        sec4(cg, r, 'RECEIVER  (band + filter + array)'); r = r+1;
        h.Blo = numL(cg, r, 'Band lo (MHz)',  2);      h.Bhi = numR(cg, r, 'Band hi (MHz)', 6); r = r+1;
        h.Filt = dropFull(cg, r, 'Filter', { ...
            'FFT brick-wall', 'FFT cosine-taper', 'FFT Gaussian', ...
            'FIR windowed-sinc', 'FIR Kaiser', ...
            'Butterworth IIR', 'Chebyshev I IIR'}, 'FFT brick-wall'); r = r+1;
        h.Ord = numL(cg, r, 'Taps/order',   101, '%d'); h.Tap = numR(cg, r, 'Taper (%)', 10);  r = r+1;
        h.Crip = numL(cg, r, 'Cheby rip dB',  1);      h.Kbeta = numR(cg, r, 'Kaiser beta', 8); r = r+1;
        h.Cest = dropFull(cg, r, 'Carrier est', ...
            {'Welch PSD', 'Periodogram (FFT)'}, 'Welch PSD'); r = r+1;
        h.Wseg = numL(cg, r, 'W. segments',   8, '%d'); h.Wov  = numR(cg, r, 'W. ovlap %', 50); r = r+1;
        h.DFm = dropFull(cg, r, 'DF mode', ...
            {'Dual-baseline', 'Single baseline'}, 'Dual-baseline'); r = r+1;
        h.Base = numL(cg, r, 'Fine base m', 0.30);     h.Bshort = numR(cg, r, 'Short base m', 0.02); r = r+1;

        sec4(cg, r, 'ARRAY CALIBRATION  (channel offsets)'); r = r+1;
        h.Phof = numL(cg, r, 'Ph off (deg)',  0);      h.Gmm = numR(cg, r, 'Gain mm (%)', 0);  r = r+1;
        h.Ucal = dropFull(cg, r, 'Use calib', {'Off', 'On'}, 'Off'); r = r+1;
        bcal = uibutton(cg, 'Text', 'CALIBRATE  (boresight reference)', ...
            'ButtonPushedFcn', @(~,~) calibrateArray());
        bcal.Layout.Row = r; bcal.Layout.Column = [1 4]; r = r+1;

        sec4(cg, r, 'TRANSMITTER  (DRFM + deception)'); r = r+1;
        h.Txm = dropFull(cg, r, 'Tx mode', ...
            {'DRFM (capture+shift)', 'Parametric regen'}, 'DRFM (capture+shift)'); r = r+1;
        h.Dly = numL(cg, r, 'RGPO dly us',    0);      h.Dop = numR(cg, r, 'Doppler kHz', 0);  r = r+1;
        h.Ntg = numL(cg, r, 'False tgts',     1, '%d'); h.Ftrf = numR(cg, r, 'Tx RF GHz', 2.1); r = r+1;
        h.Ft  = numL(cg, r, 'Tx IF MHz',      7);                                              r = r+1;

        sec4(cg, r, 'VICTIM ECCM / CLUTTER / DETECTION'); r = r+1;
        h.ECCM = dropFull(cg, r, 'ECCM mode', ...
            {'Off','PRI stagger','Frequency agility','Leading-edge tracker','All'}, 'Off'); r = r+1;
        h.VpriDev = numL(cg, r, 'Vict PRI dev %', 12); h.Vhop = numR(cg, r, 'Vict hop MHz', 1.0); r = r+1;
        h.Clutter = dropFull(cg, r, 'Ground clutter', {'Off','On'}, 'Off'); r = r+1;
        h.ClutDB = numL(cg, r, 'Clutter/skin dB', 12); h.ClutV = numR(cg, r, 'Clut spread m/s', 2); r = r+1;
        h.MTI = dropFull(cg, r, 'Clutter cancel', ...
            {'Off','2-pulse MTI','Mean cancel'}, 'Off'); r = r+1;
        h.CFAR = dropFull(cg, r, 'CFAR detection', {'Off','On'}, 'On'); r = r+1;
        h.Pfa = numL(cg, r, 'CFAR Pfa', 1e-4, '%.1e'); r = r+1;

        sec4(cg, r, 'ECM EXTENSIONS  (VGPO / cover / noise)'); r = r+1;
        h.VGPO = numL(cg, r, 'VGPO kHz/pulse', 0); h.CoverN = numR(cg, r, 'Cover pulses', 5, '%d'); r = r+1;
        h.Cover = dropFull(cg, r, 'Cover-pulse mode', {'Off','On'}, 'Off'); r = r+1;
        h.CoverLoss = numL(cg, r, 'Cover loss dB', 6); r = r+1;
        h.Jam = dropFull(cg, r, 'Noise jamming', {'Off','Spot','Barrage'}, 'Off'); r = r+1;
        h.JS = numL(cg, r, 'J/S dB', 10); h.SpotBW = numR(cg, r, 'Spot BW MHz', 0.25); r = r+1;

        sec4(cg, r, 'NOISE / REPRODUCIBILITY'); r = r+1;
        h.Seed = numL(cg, r, 'Seed',      12345, '%d');                                         r = r+1;

        % run buttons side by side
        b1 = uibutton(cg, 'Text', 'RUN PIPELINE', ...
            'ButtonPushedFcn', @(~,~) runPipeline(false));
        b1.Layout.Row = r; b1.Layout.Column = [1 2]; b1.FontWeight = 'bold';
        b2 = uibutton(cg, 'Text', 'New noise', ...
            'ButtonPushedFcn', @(~,~) runPipeline(true));
        b2.Layout.Row = r; b2.Layout.Column = [3 4]; r = r+1;

        h.results = uilabel(cg, 'Text', '', 'VerticalAlignment', 'top', ...
            'WordWrap', 'on', 'FontName', 'monospaced');
        h.results.Layout.Row = r; h.results.Layout.Column = [1 4];

        % ---------------- RIGHT: tabbed graphics ----------------
        h.tabs = uitabgroup(outer);

        t1 = uitab(h.tabs, 'Title', '1  Received');
        [h.ax_rxTime, h.ax_rxSpec] = twoAxes(t1);

        t2 = uitab(h.tabs, 'Title', '2  Filter');
        [h.ax_hmag, h.ax_himp] = twoAxes(t2);

        t3 = uitab(h.tabs, 'Title', '3  Band + Envelope');
        [h.ax_filt, h.ax_env] = twoAxes(t3);

        t4 = uitab(h.tabs, 'Title', '4  Carrier estimate');
        gg = uigridlayout(t4, [1 1]);
        h.ax_freq = uiaxes(gg);

        t5 = uitab(h.tabs, 'Title', '5  Transmit');
        [h.ax_txTime, h.ax_txSpec] = twoAxes(t5);

        t6 = uitab(h.tabs, 'Title', '6  Pipeline overview');
        og = uigridlayout(t6, [2 3]);
        h.ax_o = gobjects(1,6);
        for k = 1:6, h.ax_o(k) = uiaxes(og); end

        t7 = uitab(h.tabs, 'Title', '7  Geolocation');
        g7 = uigridlayout(t7, [1 2]);
        g7.ColumnWidth = {'1.4x', '1x'};
        h.ax_ppi  = uiaxes(g7);   % radar-scope: source position (range, azimuth)
        h.ax_chan = uiaxes(g7);   % two array channels overlaid (phase offset)

        t8 = uitab(h.tabs, 'Title', '8  Time-Frequency');
        [h.ax_sgRx, h.ax_sgTx] = twoAxes(t8);  % spectrograms: received vs Tx

        t9 = uitab(h.tabs, 'Title', '9  Victim Range-Doppler');
        g9 = uigridlayout(t9, [1 1]);
        h.ax_rd = uiaxes(g9);                  % matched-filter + slow-time FFT payoff

        t10 = uitab(h.tabs, 'Title', '10 Pulse Compression');
        [h.ax_pcWide, h.ax_pcZoom] = twoAxes(t10); % wide LFM pulse -> narrow peak

        t11 = uitab(h.tabs, 'Title', '11 Emitter Realism');
        g11 = uigridlayout(t11, [2 3]);
        h.ax_real = gobjects(1,6);
        for k = 1:6, h.ax_real(k) = uiaxes(g11); end

        t12 = uitab(h.tabs, 'Title', '12 EW Duel');
        g12 = uigridlayout(t12, [2 2]);
        h.ax_duel = gobjects(1,4);
        for k = 1:4, h.ax_duel(k) = uiaxes(g12); end

        t13 = uitab(h.tabs, 'Title', '13 Multi-Emitter ESM');
        g13 = uigridlayout(t13, [3 2]);
        g13.RowHeight = {'1x','1x',190};
        h.ax_esm = gobjects(1,4);
        for k = 1:4, h.ax_esm(k) = uiaxes(g13); end
        h.esmTable = uitable(g13, 'ColumnName', ...
            {'Track','Class','Pulses','PRI us','PW us','Freq MHz','DOA deg','Range km','Agility'}, ...
            'RowName', {});
        h.esmTable.Layout.Row = 3; h.esmTable.Layout.Column = [1 2];
    end

    function runPipeline(newNoise)
        if nargin < 1, newNoise = false; end
        if newNoise
            h.Seed.Value = h.Seed.Value + 1;   % advance to a new realization
        end
        P = readInputs();            % raw GUI -> SI parameter struct
        [P, notes] = validate(P);    % clamp / warn BEFORE computing anything

        [S.t, S.gate, S.clean, S.rx, S.rx2, S.rx3, S.emitter] = makeSignal(P, P.seed);
        R = analyzeRx(S.rx, S.rx2, S.rx3, P, S.emitter);
        S = mergeStruct(S, R);
        S.P = P; S.notes = notes;

        multiEmitterESMChain();
        transmitChain();
        victimRadarChain();
        drawAll();
        updateResults();
    end

    function P = readInputs()
        % Pull raw values from the GUI and convert GUI units (MHz, us) to SI.
        P.Fs   = h.Fs.Value  * 1e6;
        P.PRI  = h.PRI.Value * 1e-6;
        P.PW   = h.PW.Value  * 1e-6;
        P.Np   = max(1, round(h.Np.Value));
        P.SNR  = h.SNR.Value;
        P.Blo  = h.Blo.Value * 1e6;
        P.Bhi  = h.Bhi.Value * 1e6;
        P.Ft   = h.Ft.Value  * 1e6;
        P.seed = max(0, round(h.Seed.Value));

        % --- RF front-end / processing domain (new) -----------------------
        % The true carrier is at RF (GHz), but we NEVER sample at RF. A mixer
        % downconverts it to a low IF (or to complex baseband), and the ADC
        % samples there -- so Fs only has to cover the IF/baseband bandwidth,
        % which keeps the record short and the FFTs fast even for a GHz carrier.
        P.Frf  = h.Frf.Value * 1e9;            % true RF carrier (Hz), for lambda
        P.dom  = h.Dom.Value;                  % 'IF (real)' or 'Baseband IQ'
        P.iq   = strcmp(P.dom, 'Baseband IQ'); % complex baseband flag
        P.Fif  = h.Fif.Value * 1e6;            % IF centre (Hz)
        if P.iq, P.Fctr = 0; else, P.Fctr = P.Fif; end  % processing centre
        P.Flo  = P.Frf - P.Fctr;               % downconversion LO (known to Rx)
        P.Fc   = P.Fctr;                        % processing-domain carrier
        P.c    = 299792458;
        P.filt = h.Filt.Value;
        P.ord  = max(1, round(h.Ord.Value));
        P.tap  = h.Tap.Value / 100;            % percent -> fraction of BW
        P.cest = h.Cest.Value;                 % carrier estimator method
        P.wseg = max(1, round(h.Wseg.Value));  % Welch: number of segments
        P.wov  = h.Wov.Value / 100;            % Welch: overlap fraction
        P.crip = h.Crip.Value;                 % Chebyshev I passband ripple (dB)
        P.kbeta = h.Kbeta.Value;               % Kaiser window beta
        P.txmode = h.Txm.Value;                % retransmit mode (DRFM / parametric)
        % intra-pulse modulation (generator)
        P.modtype = h.Mod.Value;               % 'CW' or 'LFM chirp'
        P.bchirp  = h.Bch.Value * 1e6;         % LFM sweep bandwidth (Hz)
        % DRFM deception (transmitter)
        P.drdelay = h.Dly.Value * 1e-6;        % range gate pull-off delay (s)
        P.dop     = h.Dop.Value * 1e3;         % initial deception Doppler (Hz)
        P.ntgt    = max(1, round(h.Ntg.Value)); % number of false targets
        P.vgstep  = h.VGPO.Value * 1e3;        % VGPO increment per pulse (Hz/pulse)
        P.cover   = strcmp(h.Cover.Value,'On');
        P.coverN  = max(1,round(h.CoverN.Value));
        P.coverLoss = h.CoverLoss.Value;
        P.coverAmp = 10^(-P.coverLoss/20);
        P.jammode = h.Jam.Value;
        P.jsdb    = h.JS.Value;
        P.spotbw  = h.SpotBW.Value*1e6;

        % Victim-radar countermeasures and detector.
        P.eccm = h.ECCM.Value;
        P.eccmPRI = strcmp(P.eccm,'PRI stagger') || strcmp(P.eccm,'All');
        P.eccmFH  = strcmp(P.eccm,'Frequency agility') || strcmp(P.eccm,'All');
        P.eccmLE  = strcmp(P.eccm,'Leading-edge tracker') || strcmp(P.eccm,'All');
        P.vpriDev = h.VpriDev.Value/100;
        P.vhopstep = h.Vhop.Value*1e6;
        P.clutter = strcmp(h.Clutter.Value,'On');
        P.clutterDB = h.ClutDB.Value;
        P.clutterV = h.ClutV.Value;
        P.mtiMode = h.MTI.Value;
        P.cfar = strcmp(h.CFAR.Value,'On');
        P.pfa = h.Pfa.Value;
        % IF -> RF upconversion for transmission. The DRFM output lives at the
        % transmit IF (P.Ft, sampled/plottable); a Tx mixer with LO = FloTx
        % translates it up to the RF carrier P.Ftrf (GHz). As with the Rx side
        % we never sample at RF -- P.Ftrf is carried as a scalar and the samples
        % remain the baseband-equivalent at the IF.
        P.Ftrf  = h.Ftrf.Value * 1e9;          % transmit RF carrier target (Hz)
        P.FloTx = P.Ftrf - P.Ft;               % upconversion local oscillator (Hz)

        % --- source geometry + array (new) --------------------------------
        P.Rng  = h.Rng.Value * 1e3;            % km -> m
        P.Az   = h.Az.Value;                   % degrees (from array boresight)
        P.vrad = h.Vrad.Value;                 % m/s, positive = closing / approaching
        P.multipath = strcmp(h.MPmode.Value, 'On');
        P.mploss = h.MPloss.Value;             % reflected-path voltage loss in dB
        P.mpdelay = h.MPdly.Value * 1e-6;      % excess reflected-path delay (s)
        P.mpaz = h.MPaz.Value;                 % apparent arrival angle of reflection
        P.primode = h.PRImode.Value;
        P.pridev = h.PRIdev.Value / 100;       % fractional stagger/jitter depth
        P.fhmode = h.FHmode.Value;
        P.hopstep = h.Hstep.Value * 1e6;       % spacing between hop channels (Hz)
        % Independent multi-emitter ESM scene. The original single-emitter
        % chain remains untouched; this path builds a simultaneous three-channel
        % intercept for MUSIC, PDW deinterleaving and heuristic classification.
        P.esmOn = strcmp(h.ESMmode.Value,'On');
        P.esmN = max(2,round(h.ESMN.Value));
        P.esmScenario = h.ESMscenario.Value;
        P.musicStep = h.MusicStep.Value;
        P.pdwMiss = h.PDWMiss.Value/100;
        P.base = h.Base.Value;                 % FINE (accurate) array baseline (m)
        P.bshort = h.Bshort.Value;             % SHORT auxiliary baseline (m)
        P.dfmode = h.DFm.Value;                % 'Dual-baseline' or 'Single baseline'
        % array channel imperfections + calibration
        P.phoff  = h.Phof.Value * pi/180;      % per-channel phase-offset scale (rad)
        P.gainmm = h.Gmm.Value / 100;          % per-channel gain-mismatch scale (frac)
        P.usecal = strcmp(h.Ucal.Value, 'On'); % apply stored calibration?
    end

    function R = analyzeRx(rx1, rx2, rx3, P, M)
        % The complete receive chain (single-channel measurements use channel 1;
        % azimuth uses inter-channel phase, resolving the fine baseline's wrap
        % with the short auxiliary baseline when DF mode is 'Dual-baseline').
        Fs = P.Fs;

        % (1) band-limit ALL channels with the SELECTED filter
        R.filt  = applyFilter(rx1, P);
        filt2   = applyFilter(rx2, P);
        filt3   = applyFilter(rx3, P);

        % (2) carrier estimate (in the PROCESSING domain), refined sub-bin.
        [R.Fc_est, R.fBand, R.magBand] = estimateCarrier(R.filt, P);
        [R.fFull, R.magFull] = spectrum(R.filt, Fs);

        % Power-weighted CENTRE frequency of the occupied band (robust for a
        % wideband LFM chirp, where the peak wanders). This is what we map back
        % to RF and use for the DF wavelength.
        pb = 10.^(R.magBand/10);                  % in-band PSD (linear)
        if sum(pb) > 0
            R.Fc_ctr = sum(R.fBand .* pb) / sum(pb);
        else
            R.Fc_ctr = R.Fc_est;
        end
        % Map the processing-domain estimate back to the true RF carrier using
        % the known downconversion LO:  f_RF = f_LO + f_processing.
        R.Fc_rf = P.Flo + R.Fc_est;

        % (3) envelope -> pulse edges -> PW / PRI.
        %     Baseband is already complex, so its envelope is |x| directly; an
        %     IF signal is real, so we form the analytic signal first. z1 (the
        %     complex signal) is also what direction-finding needs.
        if P.iq
            z1 = R.filt;                  % complex baseband
        else
            z1 = analytic(R.filt);        % analytic of the real IF signal
        end
        env = abs(z1);
        w = max(1, round(0.15 * P.PW * Fs));      % smooth over ~15% of a pulse
        env = movmean(env, w);
        R.env = env;

        % Half-max threshold between the noise floor (median) and the peak. For
        % a rectangular pulse the 50% crossings approximate the true edges.
        noiseFloor = median(env);
        peak = max(env);
        R.thr = noiseFloor + 0.5 * (peak - noiseFloor);

        % Contiguous above-threshold runs are raw detections. A delayed
        % reflection can create a second run inside every PRI, so select the
        % strongest well-separated train for PW/PRI/carrier measurements while
        % retaining all runs for the envelope display.
        b = env > R.thr;
        d = diff([0 b 0]);
        riseAll = find(d == 1);
        fallAll = find(d == -1)-1;
        na = min(numel(riseAll),numel(fallAll));
        riseAll = riseAll(1:na); fallAll = fallAll(1:na);
        [rise,fall] = selectMainPulseTrain(riseAll,fallAll,env,P);
        n = min(numel(rise),numel(fall));
        R.riseAll = riseAll; R.fallAll = fallAll;
        R.rise = rise(1:n); R.fall = fall(1:n);
        R.NdetAll = na; R.Ndet = n;

        if n >= 1, R.PW_est = mean(R.fall-R.rise+1)/Fs; else, R.PW_est = NaN; end
        if n >= 2
            R.PRI_seq_est = diff(R.rise)/Fs;
            R.PRI_est = mean(R.PRI_seq_est);
            R.PRI_std = std(R.PRI_seq_est);
        else
            R.PRI_seq_est = [];
            R.PRI_est = NaN; R.PRI_std = NaN;
        end

        % Per-pulse carrier tracking, nonuniform pulse-Doppler processing, MTI,
        % per-pulse DF and the direct/multipath excess-range profile.
        A = analyzeEmitterRealism(z1, filt2, filt3, R.rise, R.fall, P, M);
        R = mergeStruct(R,A);

        % (4) LFM PULSE COMPRESSION + GEOLOCATION ------------------------
        % The envelope edge supplies a coarse TOA for every waveform. For LFM,
        % correlate the complex baseband signal with the known reference chirp.
        % The matched-filter peak is much narrower than PW (approximately 1/B),
        % so it replaces the coarse edge as the final first-pulse TOA estimate.
        if n >= 1
            R.tau_env = (R.rise(1) - 1) / Fs;      % sample 1 -> t = 0
            R.Rng_env = P.c * R.tau_env;
        else
            R.tau_env = NaN; R.Rng_env = NaN;
        end

        C = pulseCompressionRx(z1, R.Fc_ctr, P, R.rise, R.fcPulse_est);
        R.pc_active       = C.active;
        R.pc_t            = C.t;
        R.pc_mag          = C.mag;
        R.pc_magDB        = C.magDB;
        R.pc_peakIdx      = C.peakIdx;
        R.pc_peakTimes    = C.peakTimes;
        R.pc_firstPeakIdx = C.firstPeakIdx;
        R.pc_width3dB     = C.width3dB;
        R.pc_TB           = C.TB;
        R.pc_gainTheory   = C.gainTheory;
        R.pc_toaRes       = C.toaRes;
        R.pc_radarRes     = C.radarRes;
        R.pc_snrIn        = C.snrIn;
        R.pc_snrOut       = C.snrOut;
        R.pc_gainMeasured = C.gainMeasured;
        R.pc_Ncoh         = C.Ncoh;

        if C.active && isfinite(C.tauEst)
            R.tau_est = C.tauEst;
            R.Rng_est = P.c * R.tau_est;          % one-way emitter geolocation
            R.rangeMethod = 'LFM matched-filter peak';
        else
            R.tau_est = R.tau_env;
            R.Rng_est = R.Rng_env;
            R.rangeMethod = 'envelope leading edge';
        end

        % AZIMUTH by interferometry. The inter-channel phase is measured over
        % the detected-pulse samples and inverted with the RF wavelength
        % lambda = c / (f_LO + occupied-band centre).
        %
        %   FINE  baseline (element #2, x = base)   -> accurate but its phase
        %         wraps once base > lambda/2 (spatial aliasing -> wrong angle).
        %   SHORT baseline (element #3, x = bshort) -> if <= lambda/2 its phase
        %         never wraps, giving an unambiguous (but coarser) angle that is
        %         used to pick the correct wrap integer k of the fine baseline.
        if P.iq, z2 = filt2; z3 = filt3; else, z2 = analytic(filt2); z3 = analytic(filt3); end
        onmask = false(1, numel(env));
        for k = 1:n, onmask(R.rise(k):R.fall(k)) = true; end

        R.dphi_est = NaN; R.dphiS_est = NaN;
        R.Az_est = NaN; R.Az_coarse = NaN; R.Az_wrapped = false; R.kwrap = 0;
        R.Az_raw = NaN; R.cal_on = false;
        if any(onmask) && ~isnan(R.Fc_ctr) && P.base > 0
            lambda_est = P.c / (P.Flo + R.Fc_ctr);             % RF wavelength

            % measured (wrapped) phase differences, sign: angle = -dphi_geo
            dphiF = -angle(sum(z1(onmask) .* conj(z2(onmask))));  % fine
            dphiS = -angle(sum(z1(onmask) .* conj(z3(onmask))));  % short
            R.dphi_est = dphiF; R.dphiS_est = dphiS;

            % UNCALIBRATED azimuth (raw phases, still carrying any hardware
            % offset -> biased when channels are miscalibrated).
            R.Az_raw = resolveAz(dphiF, dphiS, lambda_est, P);

            % CALIBRATED: subtract the stored boresight offsets, then resolve.
            R.cal_on = P.usecal && cal.valid;
            if R.cal_on
                dphiF = dphiF - cal.dphiF;
                dphiS = dphiS - cal.dphiS;
            end
            [R.Az_est, R.Az_coarse, R.kwrap, R.Az_wrapped] = ...
                resolveAz(dphiF, dphiS, lambda_est, P);
        end
        R.filt2 = filt2; R.filt3 = filt3;
    end

    function A = analyzeEmitterRealism(z1,filt2,filt3,rise,fall,P,M)
        % Measurements behind the Emitter Realism tab. Carrier is tracked on
        % each pulse after optional LFM dechirping. Slow-time Doppler uses a
        % nonuniform DFT, so staggered/jittered PRI samples need not be forced
        % onto a fictitious uniform grid.
        A = struct('fcPulse_est',[],'fcPulse_true',[],'hop_est',[], ...
            'fc_center_est',NaN,'carrier_rmse',NaN,'pri_true',M.priTrue, ...
            'vel_axis',[],'doppler_db',[],'fd_rx_est',NaN,'vel_rx_est',NaN, ...
            'slow_amp',[],'mti_amp',[],'mti_gain',NaN,'az_perPulse',[], ...
            'az_spread',NaN,'range_excess',[],'range_profile_db',[]);
        if isempty(rise), return; end
        if P.iq, z2=filt2; z3=filt3; else, z2=analytic(filt2); z3=analytic(filt3); end
        Fs=P.Fs; L=max(8,round(P.PW*Fs)); np=numel(rise);
        fhat=NaN(1,np); coherent=NaN(1,np); azp=NaN(1,np);
        fcTrue=NaN(1,np);
        mapIdx=zeros(1,np);
        for k=1:np
            te=(rise(k)-1)/Fs;
            [~,mapIdx(k)]=min(abs(M.arrivalTimes-te));
            fcTrue(k)=M.fcProcTrue(mapIdx(k));
        end

        % First pass: dechirped per-pulse FFT carrier estimates.
        Nfft=2^nextpow2(max(2048,8*L));
        fax=((0:Nfft-1)-Nfft/2)/Nfft*Fs;
        for k=1:np
            ia=rise(k); ib=min(numel(z1),ia+L-1); idx=ia:ib;
            tt=(idx-1)/Fs; tl=tt-((rise(k)-1)/Fs+P.PW/2);
            seg=z1(idx);
            if strcmp(P.modtype,'LFM chirp') && P.bchirp>0
                seg=seg.*exp(-1i*pi*(P.bchirp/P.PW)*tl.^2);
            end
            w=0.5-0.5*cos(2*pi*(0:numel(seg)-1)/max(1,numel(seg)-1));
            X=abs(fftshift(fft(seg.*w,Nfft)));
            if P.iq
                mask=abs(fax-P.Fctr)<=max(P.Fs/3,3*P.hopstep+abs(P.fdRx)+P.bchirp);
            else
                mask=fax>=0 & fax<=P.Fs/2;
            end
            xm=X; xm(~mask)=0; [~,im]=max(xm);
            delta=0;
            if im>1 && im<Nfft
                a=log(X(im-1)+eps); b=log(X(im)+eps); c=log(X(im+1)+eps);
                den=a-2*b+c;
                if abs(den)>eps, delta=max(-0.5,min(0.5,0.5*(a-c)/den)); end
            end
            fhat(k)=fax(im)+delta*Fs/Nfft;
        end

        % Quantise hop channels about the nominal processing centre; the common
        % residual after removing those channels is receive Doppler plus estimator
        % bias, so hopping is removed without erasing slow-time phase progression.
        if P.hopstep>0
            lev=round((fhat-P.Fctr)/P.hopstep);
            if strcmp(P.fhmode,'3-tone hop'), lev=max(-1,min(1,lev));
            elseif strcmp(P.fhmode,'Random hop'), lev=max(-2,min(2,lev)); end
            hopEst=lev*P.hopstep;
        else
            hopEst=zeros(size(fhat));
        end
        common=median(fhat-hopEst);

        for k=1:np
            ia=rise(k); ib=min(numel(z1),ia+L-1); idx=ia:ib;
            tt=(idx-1)/Fs; tl=tt-((rise(k)-1)/Fs+P.PW/2);
            ref=ones(size(tt));
            if strcmp(P.modtype,'LFM chirp') && P.bchirp>0
                ref=exp(1i*pi*(P.bchirp/P.PW)*tl.^2);
            end
            dehop=exp(-1i*2*pi*(P.Fctr+hopEst(k))*tt);
            coherent(k)=sum(z1(idx).*dehop.*conj(ref))/numel(idx);

            dF=-angle(sum(z1(idx).*conj(z2(idx))));
            dS=-angle(sum(z1(idx).*conj(z3(idx))));
            if P.usecal && cal.valid
                dF=dF-cal.dphiF; dS=dS-cal.dphiS;
            end
            lambda=P.c/max(P.Flo+fhat(k),eps);
            azp(k)=resolveAz(dF,dS,lambda,P);
        end

        slowT=(rise-1)/Fs;
        if np>=2
            dt=diff(slowT); fmax=0.5/max(min(dt),eps);
            fgrid=linspace(-fmax,fmax,1025);
            wp=0.5-0.5*cos(2*pi*(0:np-1)/max(1,np-1));
            E=exp(-1i*2*pi*(fgrid(:)*(slowT-slowT(1))));
            pd=abs(E*(coherent.*wp).').^2;
            pdDB=10*log10(pd/(max(pd)+eps)+eps); pdDB=max(pdDB,-60);
            [~,ip]=max(pd); fdHat=fgrid(ip);
            velAxis=fgrid*P.c/max(P.Frf,eps);
            velHat=fdHat*P.c/max(P.Frf,eps);
            mti=[0 diff(coherent)];
            A.mti_gain=20*log10((sqrt(mean(abs(mti).^2))+eps)/(sqrt(mean(abs(coherent).^2))+eps));
        else
            fgrid=[]; pdDB=[]; fdHat=NaN; velAxis=[]; velHat=NaN; mti=0;
        end

        % First-pulse excess-path range profile. LFM gives a compressed direct
        % peak and a second peak at c*extra_delay; CW truthfully remains broad.
        pre=round(0.15*P.PW*Fs);
        post=round((P.PW+P.mpdelay+0.5*P.PW)*Fs);
        ia=max(1,rise(1)-pre); ib=min(numel(z1),rise(1)+post); idx=ia:ib;
        tt=(idx-1)/Fs;
        if strcmp(P.modtype,'LFM chirp') && P.bchirp>0
            tl=(0:L-1)/Fs-P.PW/2;
            ref=exp(1i*pi*(P.bchirp/P.PW)*tl.^2);
            bb=z1(idx).*exp(-1i*2*pi*fhat(1)*tt);
            mf=conv(bb,conj(fliplr(ref)),'full'); rp=abs(mf);
            tprof=((ia-1)+(1:numel(mf))-L)/Fs;
        else
            rp=abs(z1(idx)); tprof=(idx-1)/Fs;
        end
        rex=P.c*(tprof-(rise(1)-1)/Fs);
        rpDB=20*log10(rp/(max(rp)+eps)+eps); rpDB=max(rpDB,-70);

        A.fcPulse_est=fhat; A.fcPulse_true=fcTrue; A.hop_est=hopEst;
        A.fc_center_est=common; A.carrier_rmse=sqrt(mean((fhat-fcTrue).^2));
        A.vel_axis=velAxis; A.doppler_db=pdDB; A.fd_rx_est=fdHat;
        A.vel_rx_est=velHat; A.slow_amp=coherent; A.mti_amp=mti;
        A.az_perPulse=azp; A.az_spread=std(azp(isfinite(azp)));
        A.range_excess=rex; A.range_profile_db=rpDB;
    end

    function calibrateArray()
        % Measure the array's channel phase (and gain) offsets from a BORESIGHT
        % reference. At azimuth 0 the geometric phase difference is zero, so any
        % measured inter-channel phase is pure hardware offset -- we store it and
        % subtract it during direction finding. (Real systems inject a common
        % reference tone or use a boresight source; same idea.)
        P = validate(readInputs());
        Pc=P; Pc.Az=0; Pc.SNR=max(P.SNR,40); Pc.modtype='CW'; Pc.bchirp=0;
        Pc.vrad=0; Pc.fdRx=0; Pc.multipath=false; Pc.mpamp=0;
        Pc.primode='Constant'; Pc.pridev=0; Pc.fhmode='Fixed'; Pc.hopstep=0;
        [~, ~, ~, r1, r2, r3] = makeSignal(Pc, P.seed);
        f1 = applyFilter(r1, Pc); f2 = applyFilter(r2, Pc); f3 = applyFilter(r3, Pc);
        if Pc.iq, z1 = f1; z2 = f2; z3 = f3;
        else,     z1 = analytic(f1); z2 = analytic(f2); z3 = analytic(f3); end
        [rs, fs] = detectPulses(f1, Pc.PW, Pc.Fs);
        on = false(1, numel(z1));
        for k = 1:numel(rs), on(rs(k):fs(k)) = true; end
        if ~any(on), on(:) = true; end
        cal.dphiF = -angle(sum(z1(on) .* conj(z2(on))));   % = hardware offset (fine)
        cal.dphiS = -angle(sum(z1(on) .* conj(z3(on))));   % = hardware offset (short)
        a1 = mean(abs(z1(on)));
        cal.gF = mean(abs(z2(on))) / max(a1, eps);
        cal.gS = mean(abs(z3(on))) / max(a1, eps);
        cal.valid = true;
        % turn calibration on and refresh so the user sees the corrected result
        h.Ucal.Value = 'On';
        runPipeline(false);
    end

% ---------------------------------------------------------- Tx CHAIN
    function transmitChain()
        % Two retransmit modes:
        %
        %   'DRFM (capture+shift)'  -- Digital RF Memory. The receiver STORES the
        %       actual captured samples (here the band-limited channel-1 record)
        %       and replays them, frequency-translated to the target frequency.
        %       Because it replays the real waveform, it PRESERVES the pulse
        %       envelope, exact timing/delay, and any intra-pulse modulation
        %       (chirp, phase codes) -- the whole reason DRFM exists. The
        %       translation amount uses the ESTIMATED carrier: df = Ft - fc_hat.
        %
        %   'Parametric regen'      -- the old behaviour: synthesise a CLEAN
        %       pulse train from the MEASURED PW/PRI and modulate a fresh carrier
        %       at Ft. Discards intra-pulse structure and the capture noise.
        P = S.P; Fs = P.Fs;

        if strcmp(P.txmode, 'Parametric regen')
            if ~strcmp(P.primode,'Constant') || ~strcmp(P.fhmode,'Fixed')
                S.notes{end+1} = ['Parametric regeneration uses mean PW/PRI and a single carrier; ' ...
                    'it does not preserve emitter PRI/frequency agility. Select DRFM to replay it.'];
            end
            PW  = S.PW_est;  if isnan(PW),  PW  = P.PW;  end
            PRI = S.PRI_est; if isnan(PRI), PRI = P.PRI; end
            if S.Ndet >= 1, Np = S.Ndet; else, Np = P.Np; end
            Ntot = round(Np * PRI * Fs) + round(0.5*PRI*Fs);
            t = (0:Ntot-1) / Fs;
            gate = zeros(1, Ntot);
            for k = 0:Np-1
                t0 = k * PRI;
                gate(t >= t0 & t < t0 + PW) = 1;
            end
            % modulate at Ft in the processing domain: real cos for IF, complex
            % exponential for baseband.
            if P.iq
                S.tx = gate .* exp(1i*2*pi*P.Ft*t);
            else
                S.tx = gate .* cos(2*pi*P.Ft*t);
            end
            S.tx_t     = t;
            S.tx_shift = NaN;
            S.ft_delays = [];
            S.tx_ghostR = []; S.tx_ghostW = [];
            S.drfm_base = [];
            S.cover_delays = [];
        else
            % ---- DRFM: replay the stored capture, with deception ----------
            % df = (Ft - fc_hat) + Doppler. For IF the stored capture is real, so
            % we translate its analytic signal and keep the real part; for
            % baseband it is already complex and we translate it directly. Either
            % way the envelope, timing and intra-pulse modulation are preserved.
            % Deception adds a range gate pull-off (delay) and several delayed
            % copies as MULTIPLE FALSE TARGETS.
            x = S.filt;                      % stored capture, ch1
            if isfield(S,'fc_center_est') && isfinite(S.fc_center_est)
                fc_hat = S.fc_center_est;
            else
                fc_hat = S.Fc_est;
            end
            if isnan(fc_hat), fc_hat = P.Fc; end
            df = (P.Ft - fc_hat) + P.dop;    % translation incl. Doppler
            n  = 0:numel(x)-1;
            if P.iq
                base = x .* exp(1i*2*pi*df*n/Fs);        % complex baseband
            else
                base = real(analytic(x) .* exp(1i*2*pi*df*n/Fs));
            end
            base = applyVGPO(base,S.rise,S.fall,P);
            S.drfm_base = base;                         % actual replay train before delays

            PW = S.PW_est; if isnan(PW), PW = P.PW; end
            N   = max(1, round(P.ntgt));
            % Target spacing = the RGPO delay. If it is smaller than a pulse
            % width the false targets OVERLAP and their pulses interfere -- the
            % detector below then merges them, which is shown truthfully on the
            % PPI. (With no RGPO, fall back to a clean 2*PW separation.)
            if P.drdelay > 0, step = P.drdelay; else, step = 2*PW; end
            dmax = round((P.drdelay + (N-1)*step) * Fs);
            y = zeros(1, numel(base) + dmax);
            if ~isreal(base), y = complex(y); end
            delays = zeros(1, N);
            for m = 0:N-1
                delays(m+1) = P.drdelay + m*step;
                d = round(delays(m+1) * Fs);
                seg = d + (1:numel(base));
                y(seg) = y(seg) + base;      % equal-amplitude false targets
            end
            [y,S.cover_delays] = addCoverPulses(y,base,P,P.seed+701);
            y = addNoiseJammer(y,P,P.seed+1701);
            S.tx        = y;
            S.tx_t      = (0:numel(y)-1) / Fs;
            S.tx_shift  = df;
            S.ft_delays = delays;            % nominal (theoretical) delays

            % REAL false-target positions: run the SAME pulse detector on the
            % actual transmitted waveform. If two false targets fall within a
            % pulse width of each other their echoes MERGE (interfere) and the
            % detector returns fewer / wider pulses -- so the PPI then shows the
            % deception exactly as it comes out, overlaps and all, instead of the
            % ideal delay math.
            PRIx = S.PRI_est; if isnan(PRIx), PRIx = P.PRI; end
            [rT, fT] = detectPulses(real(S.tx), PW, Fs);
            if ~isempty(rT)
                tRise = (rT - 1) / Fs;                 % onset times (s)
                win   = tRise(1) + PRIx - 0.25*PW;     % first range window
                keep  = tRise <= win;
                S.tx_ghostR   = P.c * tRise(keep);     % detected apparent range (m)
                S.tx_ghostW   = (fT(keep) - rT(keep) + 1) / Fs;  % detected width (s)
            else
                S.tx_ghostR = []; S.tx_ghostW = [];
            end

            % Flag if the translated band would run past the sampled edge.
            halfBW = 0.5*(P.Bhi - P.Blo);
            if P.iq
                if abs(P.dop) + halfBW > Fs/2
                    S.notes{end+1} = 'DRFM: shifted baseband exceeds +/-Fs/2 (aliasing).';
                end
            else
                if P.Ft + P.dop + halfBW > Fs/2 || P.Ft + P.dop - halfBW < 0
                    S.notes{end+1} = 'DRFM: translated band exceeds [0, Nyquist] (aliasing).';
                end
            end
        end

        [S.tx_f, S.tx_mag] = spectrum(S.tx, Fs);

        % ---- IF -> RF UPCONVERSION -------------------------------------
        % A transmit mixer with local oscillator FloTx translates the sampled
        % IF signal (centred at Ft, +Doppler) up to the RF carrier. We do not
        % synthesise the GHz waveform (that would need Fs > 2*f_RF); instead the
        % sampled S.tx IS the baseband-equivalent of the RF output, and the true
        % transmitted RF carrier is recorded as a scalar:
        %     f_RF_tx = (Ft + Doppler) + FloTx = Ftrf + Doppler.
        S.tx_lo = P.FloTx;                        % upconversion LO (Hz)
        S.tx_rf = P.Ftrf + P.dop;                 % transmitted RF carrier (Hz)
    end


    function y = applyVGPO(x,rise,fall,P)
        % Pulse-to-pulse Doppler walk. The extra frequency on pulse k is
        % (k-1)*vgstep; phase is integrated between pulses so the slow-time
        % sequence is coherent rather than independently rephased.
        if abs(P.vgstep) < eps || isempty(rise)
            y=x; return
        end
        wasReal=isreal(x);
        if wasReal, z=analytic(x); else, z=x; end
        np=min([numel(rise),numel(fall),P.Np]);
        phi=0;
        for k=1:np
            a=max(1,rise(k)); b=min(numel(z),fall(k));
            if b<a, continue, end
            tl=(0:b-a)/P.Fs;
            fd=(k-1)*P.vgstep;
            z(a:b)=z(a:b).*exp(1i*(phi+2*pi*fd*tl));
            if k<np
                if isfield(S,'emitter') && numel(S.emitter.arrivalTimes)>=k+1
                    dt=S.emitter.arrivalTimes(k+1)-S.emitter.arrivalTimes(k);
                else
                    dt=P.PRI;
                end
                phi=phi+2*pi*fd*dt;
            end
        end
        if wasReal, y=real(z); else, y=z; end
    end

    function victimRadarChain()
        % Full victim-versus-jammer model. The raw map contains skin, DRFM,
        % clutter and noise-jamming returns. The defended map then applies the
        % selected PRI/frequency ECCM, leading-edge tracking, MTI and CA-CFAR.
        P=S.P; Fs=P.Fs; PRI=P.PRI; PW=P.PW;
        Np=min(P.Np,32);
        if Np<P.Np
            S.notes{end+1}=sprintf('Victim radar CPI uses the first %d of %d pulses.',Np,P.Np);
        end

        L=max(2,round(PW*Fs)); tp=(0:L-1)/Fs;
        if strcmp(P.modtype,'LFM chirp') && P.bchirp>0
            K=P.bchirp/PW; tl=tp-PW/2;
            ref=exp(1i*pi*K*tl.^2);
            S.victim_rangeRes=P.c/(2*P.bchirp);
        else
            ref=ones(1,L);
            S.victim_rangeRes=P.c*PW/2;
        end
        ref=ref/sqrt(sum(abs(ref).^2)+eps);
        hmf=conj(fliplr(ref));

        % Victim timing/frequency ECCM. Stagger is deliberately independent of
        % the intercepted emitter schedule: an unsophisticated repeater that uses
        % captured timing/tone therefore loses range/frequency coherence.
        priSeq=PRI*ones(1,max(0,Np-1));
        if P.eccmPRI && ~isempty(priSeq)
            pat=repmat([-1 1],1,ceil(numel(priSeq)/2));
            priSeq=PRI*(1+P.vpriDev*pat(1:numel(priSeq)));
        end
        tSlow=[0 cumsum(priSeq)];
        if P.eccmFH
            hopLev=mod(0:Np-1,3)-1;
            victimHop=hopLev*P.vhopstep;
        else
            victimHop=zeros(1,Np);
        end

        isDRFM=strcmp(P.txmode,'DRFM (capture+shift)');
        deceptionActive=isDRFM && (P.drdelay>0 || P.ntgt>1 || abs(P.dop)>0 || ...
            abs(P.vgstep)>0 || P.cover || ~strcmp(P.jammode,'Off'));
        ghostDelay=[];
        if deceptionActive && isfield(S,'ft_delays') && ~isempty(S.ft_delays)
            ghostDelay=S.ft_delays(:).';
        elseif deceptionActive
            ghostDelay=0;
        end
        Ng=numel(ghostDelay);

        replay=extractVictimReplayPulse(ref,L,P);
        fdPulse=P.dop+(0:Np-1)*P.vgstep;
        phiPulse=zeros(1,Np);
        for k=2:Np
            phiPulse(k)=phiPulse(k-1)+2*pi*fdPulse(k-1)*priSeq(k-1);
        end

        % Range ambiguity under PRI stagger: delays beyond a PRI wrap to a
        % different apparent range on alternating PRIs and lose coherent support.
        priForPulse=PRI*ones(1,Np);
        if ~isempty(priSeq), priForPulse=[priSeq priSeq(end)]; end
        ghostEff=zeros(Ng,Np);
        for m=1:Ng
            for k=1:Np
                if P.eccmPRI && ghostDelay(m)>=priForPulse(k)
                    ghostEff(m,k)=mod(ghostDelay(m),priForPulse(k));
                else
                    ghostEff(m,k)=ghostDelay(m);
                end
            end
        end

        coverDelays=[];
        if P.cover
            span=min(0.8*PRI,max([4*PW,P.drdelay+2*PW,2*PW]));
            coverDelays=linspace(0.15*PW,span,P.coverN);
        end
        maxEff=max([0,ghostEff(:).',coverDelays]);
        guard=max(round(0.75*L),8); tail=max(round(0.75*L),16);
        Nfast=guard+round(maxEff*Fs)+L+tail;
        Nrc=Nfast+L-1;

        skinAmp=0.42;
        victimSNR=min(35,max(8,P.SNR));
        noisePow=(1/max(L,1))/10^(victimSNR/10);

        % Streaming energy pass limits display memory for large pull-offs.
        score=zeros(Nrc,1); rs1=RandStream('mt19937ar','Seed',P.seed+4242);
        for k=1:Np
            xpri=makeVictimPRIDuel(k,replay,ref,guard,Nfast,Fs,skinAmp, ...
                ghostEff,fdPulse,phiPulse,victimHop,coverDelays,noisePow,rs1,P);
            rck=conv(xpri,hmf,'full'); score=score+abs(rck).^2;
        end
        maxRangeCells=6000;
        if Nrc>maxRangeCells
            stride=ceil(Nrc/maxRangeCells); keep=zeros(1,ceil(Nrc/stride)); q=0;
            for a=1:stride:Nrc
                q=q+1; b=min(Nrc,a+stride-1); [~,jj]=max(score(a:b)); keep(q)=a+jj-1;
            end
            keep=keep(1:q);
            S.notes{end+1}=sprintf('Victim display pooled %d samples into %d range cells.',Nrc,numel(keep));
        else
            keep=1:Nrc;
        end

        rcMap=complex(zeros(numel(keep),Np));
        rs2=RandStream('mt19937ar','Seed',P.seed+4242);
        for k=1:Np
            xpri=makeVictimPRIDuel(k,replay,ref,guard,Nfast,Fs,skinAmp, ...
                ghostEff,fdPulse,phiPulse,victimHop,coverDelays,noisePow,rs2,P);
            rck=conv(xpri,hmf,'full'); rcMap(:,k)=rck(keep);
        end

        skinPeak=guard+L;
        rangeAxis=P.Rng+P.c*(keep-skinPeak)/Fs/2;
        valid=rangeAxis>=0; rangeAxis=rangeAxis(valid); rcMap=rcMap(valid,:);

        % Track the pulse-compressed target returns before clutter is injected.
        % This models a narrow range-tracking gate/leading-edge discriminator.
        rcTrack=rcMap;
        rcRaw=rcMap;

        % Range-compressed stationary ground clutter creates the zero-Doppler
        % ridge. A small random velocity spread broadens it realistically.
        if P.clutter
            rcRaw=addVictimClutter(rcRaw,rangeAxis,tSlow,skinAmp,P);
        end

        % PRI-stagger ECCM compares the two stagger groups in range. A genuine
        % return remains in the same cell, whereas a range-ambiguous repeater
        % wraps to different cells and receives a low consistency weight.
        rcDef=rcRaw; S.victim_priCellReject=0;
        if P.eccmPRI
            [priWeight,S.victim_priCellReject]=priConsistencyWeights(rcTrack);
            rcDef=rcDef.*priWeight;
        end

        % A conventional peak tracker follows the strongest RGPO return; the
        % leading-edge tracker selects the first statistically significant rise
        % inside the predicted gate.
        [trackPeak,trackLead]=victimRangeTrackers(rcTrack,rangeAxis,P,S.victim_rangeRes,maxEff);
        if P.eccmLE, trackUsed=trackLead; else, trackUsed=trackPeak; end
        S.victim_trackPeak=trackPeak; S.victim_trackLead=trackLead;
        S.victim_trackUsed=trackUsed; S.victim_tSlow=tSlow;

        % Raw cluttered map, then selected clutter canceller.
        Ndop=max(128,2^nextpow2(8*Np));
        [RDraw,fdAxis]=slowTimeDFT(rcRaw,tSlow,P.Ftrf,PRI,Ndop);
        switch P.mtiMode
            case '2-pulse MTI'
                rcProc=diff(rcDef,1,2); tProc=tSlow(2:end);
            case 'Mean cancel'
                rcProc=rcDef-mean(rcDef,2); tProc=tSlow;
            otherwise
                rcProc=rcDef; tProc=tSlow;
        end
        [RDproc,fdAxisProc]=slowTimeDFT(rcProc,tProc,P.Ftrf,PRI,Ndop);
        fdAxis=fdAxisProc;

        rawMag=abs(RDraw); procMag=abs(RDproc);
        normRef=max([rawMag(:);procMag(:);eps]);
        S.victim_rdRawDB=20*log10(rawMag/normRef+eps);
        S.victim_rdRawDB=max(S.victim_rdRawDB,-60);
        S.victim_rdDB=20*log10(procMag/normRef+eps);
        S.victim_rdDB=max(S.victim_rdDB,-60);
        S.victim_range=rangeAxis;
        lambda=P.c/max(P.Ftrf,eps);
        S.victim_velocity=fdAxis*lambda/2;
        [~,iz]=min(abs(S.victim_velocity));
        S.victim_clutterSuppression=20*log10((max(rawMag(:,iz))+eps)/(max(procMag(:,iz))+eps));
        S.victim_velRes=lambda/(2*max(1,numel(tProc))*mean([priSeq PRI]));
        S.victim_velUnamb=lambda/(4*PRI);
        S.victim_prf=1/PRI; S.victim_snr=victimSNR; S.victim_sampleRate=Fs;

        % Expected trajectories and ECCM effects.
        S.victim_ghostRangePulse=P.Rng+P.c*ghostEff/2;
        S.victim_ghostVelPulse=repmat(fdPulse*lambda/2,Ng,1);
        if Ng>0
            S.victim_ghostRange=median(S.victim_ghostRangePulse,2).';
            S.victim_ghostVel=median(S.victim_ghostVelPulse,2).';
        else
            S.victim_ghostRange=[]; S.victim_ghostVel=[];
        end
        S.victim_fdPulse=fdPulse; S.victim_vgpoVel=fdPulse*lambda/2;
        S.victim_coverRange=P.Rng+P.c*coverDelays/2;
        S.victim_hop=victimHop; S.victim_priSeq=priSeq;
        S.victim_priReject=0;
        if Ng>0 && P.eccmPRI
            spread=std(S.victim_ghostRangePulse,0,2);
            S.victim_priReject=sum(spread>2*S.victim_rangeRes);
        end

        % 2-D cell-averaging CFAR on the defended map. Local-maximum pruning
        % turns each extended target blob into a compact detection list.
        S.victim_cfarMask=false(size(procMag)); S.victim_cfarThreshold=zeros(size(procMag));
        S.victim_detRange=[]; S.victim_detVel=[]; S.victim_detPower=[];
        if P.cfar && ~isempty(procMag)
            [mask,thr]=caCfar2D(procMag.^2,P.pfa);
            mask=mask & localMax2(procMag);
            [rr,cc]=find(mask);
            if numel(rr)>120
                vals=procMag(sub2ind(size(procMag),rr,cc));
                [~,ord]=sort(vals,'descend'); ord=ord(1:120); rr=rr(ord); cc=cc(ord);
            end
            S.victim_cfarMask=mask; S.victim_cfarThreshold=thr;
            S.victim_detRange=rangeAxis(rr(:)).'; S.victim_detRange=S.victim_detRange(:);
            S.victim_detVel=S.victim_velocity(cc(:)).'; S.victim_detVel=S.victim_detVel(:);
            S.victim_detPower=procMag(sub2ind(size(procMag),rr,cc));
        end
        [S.victim_skinDetected,S.victim_ghostDetected]=classifyVictimDetections(P);

        S.victim_dopplerAliased=any(abs(fdPulse)>(0.5/PRI));
        if S.victim_dopplerAliased
            S.notes{end+1}='Victim radar: part of the RGPO/VGPO Doppler trajectory aliases at the nominal PRF.';
        end
        if strcmp(P.modtype,'CW')
            S.notes{end+1}=sprintf('Victim radar: CW range resolution is broad (about %.2f km).',S.victim_rangeRes/1e3);
        end
    end

    function replay=extractVictimReplayPulse(ref,L,P)
        % Extract the actual filtered/noisy DRFM pulse shape, then remove its
        % nominal IF and fixed Doppler so victimRadarChain can impose arbitrary
        % RGPO/VGPO and ECCM schedules coherently.
        replay=ref(:).';
        if ~isfield(S,'drfm_base') || isempty(S.drfm_base) || isempty(S.rise), return, end
        z=S.drfm_base; if isreal(z), z=analytic(z); end
        a=max(1,S.rise(1)); b=min(numel(z),a+L-1);
        q=z(a:b); if numel(q)<L, q(end+1:L)=0; end
        n=a-1+(0:L-1); q=q.*exp(-1i*2*pi*P.Ft*n/P.Fs);
        q=q.*exp(-1i*2*pi*P.dop*(0:L-1)/P.Fs);
        if any(abs(q)>0), replay=q/sqrt(sum(abs(q).^2)+eps); end
    end

    function [skinHit,ghostHit]=classifyVictimDetections(P)
        skinHit=false; ghostHit=0;
        if isempty(S.victim_detRange), return, end
        dr=max(2*S.victim_rangeRes,50); dv=max(2*S.victim_velRes,1);
        skinHit=any(abs(S.victim_detRange-P.Rng)<=dr & abs(S.victim_detVel)<=dv);
        for m=1:numel(S.victim_ghostRange)
            rr=S.victim_ghostRangePulse(m,:); vv=S.victim_ghostVelPulse(m,:);
            hit=false;
            for k=1:numel(rr)
                if any(abs(S.victim_detRange-rr(k))<=dr & abs(S.victim_detVel-vv(k))<=dv)
                    hit=true; break
                end
            end
            ghostHit=ghostHit+hit;
        end
    end


% =====================================================================
%                MULTI-EMITTER ESM / MUSIC / DEINTERLEAVING
% =====================================================================
    function multiEmitterESMChain()
        % Build an independent simultaneous-emitter intercept. The main app's
        % single-emitter signal continues to drive the DRFM/EW-duel path; this
        % second scene demonstrates the front end of a practical ESM receiver:
        % three-channel acquisition -> PDWs -> deinterleaving -> classification.
        P = S.P;
        S.esm = emptyESMResult();
        if ~P.esmOn, return; end

        E = buildESMProfiles(P);
        [X,t,pulseTruth,noiseStd] = simulateESMScene(P,E,P.seed+19000);
        pdw = extractESMPDWs(P,X,t,pulseTruth,noiseStd,P.seed+20000);
        [labels,tracks] = deinterleaveESMPDWs(P,pdw,E);
        music = runESMMUSIC(P,X,t,pdw,labels,tracks);

        % Replace rough beamformer DOA with the nearest group MUSIC peak when
        % available, then refresh classification/table statistics.
        if ~isempty(tracks) && ~isempty(music.trackDoa)
            for k=1:min(numel(tracks),numel(music.trackDoa))
                if isfinite(music.trackDoa(k)), tracks(k).doa = music.trackDoa(k); end
            end
        end

        [sgT,sgF,sgDB] = esmSTFT(X(1,:),P.Fs,P.iq);
        R = emptyESMResult();
        R.active = true; R.emitters = E; R.X = X; R.t = t;
        R.truth = pulseTruth; R.pdw = pdw; R.labels = labels;
        R.tracks = tracks; R.music = music;
        R.sgT = sgT; R.sgF = sgF; R.sgDB = sgDB;
        R.detectedPDW = numel(pdw); R.truePulses = numel(pulseTruth);
        if isempty(tracks), R.meanPurity=NaN;
        else, R.meanPurity=mean([tracks.purity]); end
        S.esm = R;
    end

    function drawMultiEmitterESM(P)
        ax=h.ax_esm;
        for k=1:4, cla(ax(k)); end
        if ~isfield(S,'esm') || ~S.esm.active
            for k=1:4
                text(ax(k),0.5,0.5,'Turn ESM scene On in the controls', ...
                    'Units','normalized','HorizontalAlignment','center'); axis(ax(k),'off');
            end
            h.esmTable.Data={'--','ESM disabled',0,NaN,NaN,NaN,NaN,NaN,'--'};
            return
        end
        R=S.esm; MHz=1e-6; us=1e6;
        for k=1:4, axis(ax(k),'on'); end
        imagesc(ax(1),R.sgT*us,R.sgF*MHz,R.sgDB); set(ax(1),'YDir','normal');
        try, colormap(ax(1),'turbo'); catch, colormap(ax(1),'jet'); end
        title(ax(1),sprintf('Composite intercept: %d emitters / %d PDWs',numel(R.emitters),numel(R.pdw)));
        xlabel(ax(1),'Time (\mus)'); ylabel(ax(1),'Processing frequency (MHz)');
        if ~P.iq, ylim(ax(1),[max(0,P.Blo*MHz-0.5) min(P.Fs/2*MHz,P.Bhi*MHz+0.5)]); end

        hold(ax(2),'on'); leg={}; lh=gobjects(1,0);
        for g=1:numel(R.music.groups)
            z=R.music.groups(g); lh(end+1)=plot(ax(2),R.music.grid,z.pdb,'LineWidth',1.2); %#ok<AGROW>
            leg{end+1}=sprintf('%.3f MHz, tracks %s',z.freq*MHz,mat2str(z.members)); %#ok<AGROW>
            for q=1:numel(z.peaks), xline(ax(2),z.peaks(q),'--','HandleVisibility','off'); end
        end
        for e=1:numel(R.emitters), xline(ax(2),R.emitters(e).az,':','HandleVisibility','off'); end
        hold(ax(2),'off'); grid(ax(2),'on'); ylim(ax(2),[-55 3]); xlim(ax(2),[-90 90]);
        title(ax(2),'Per-frequency MUSIC pseudospectra (dotted = true DOA)');
        xlabel(ax(2),'Azimuth (deg)'); ylabel(ax(2),'Normalized spectrum (dB)');
        if ~isempty(leg), legend(ax(2),lh,leg,'Location','best'); end

        if ~isempty(R.pdw)
            lab=R.labels(:); scatter(ax(3),[R.pdw.toa]*us,[R.pdw.freq]*MHz,28,lab,'filled');
            title(ax(3),'Merged PDW stream after deinterleaving');
            xlabel(ax(3),'TOA (\mus)'); ylabel(ax(3),'Measured frequency (MHz)'); grid(ax(3),'on');
        end

        T=R.tracks; if ~isempty(T)
            scatter(ax(4),[T.PRI]*us,[T.pw]*us,72,1:numel(T),'filled'); hold(ax(4),'on');
            for k=1:numel(T), text(ax(4),T(k).PRI*us,T(k).pw*us,['  T' num2str(k) ': ' T(k).class]); end
            hold(ax(4),'off'); grid(ax(4),'on');
            title(ax(4),'Track PRI/PW signatures and heuristic classification');
            xlabel(ax(4),'Estimated PRI (\mus)'); ylabel(ax(4),'Estimated PW (\mus)');
            data=cell(numel(T),9);
            for k=1:numel(T)
                data(k,:)={sprintf('T%d',k),T(k).class,T(k).count,T(k).PRI*us, ...
                    T(k).pw*us,T(k).freqMed*MHz,T(k).doa,T(k).range/1e3,T(k).agility};
            end
            h.esmTable.Data=data;
        else
            h.esmTable.Data={'--','No tracks',0,NaN,NaN,NaN,NaN,NaN,'--'};
        end
    end

    function txt=esmResultsText()
        if ~isfield(S,'esm') || ~S.esm.active
            txt=sprintf('\n\nMULTI-EMITTER ESM\n-------------------------\nStatus       : OFF'); return
        end
        R=S.esm; P=S.P;
        txt=sprintf(['\n\nMULTI-EMITTER ESM\n-------------------------\n' ...
            'Scenario     : %s\nEmitters     : %d\nTrue pulses  : %d\nMeasured PDWs: %d\n' ...
            'Tracks       : %d\nMean purity  : %.1f %%\nMUSIC groups : %d\n' ...
            'Array limit  : max 2 co-channel sources with 3 sensors\n'], ...
            P.esmScenario,numel(R.emitters),R.truePulses,R.detectedPDW, ...
            numel(R.tracks),100*R.meanPurity,numel(R.music.groups));
        for k=1:numel(R.tracks)
            T=R.tracks(k);
            txt=[txt sprintf('T%d %-20s PRI %7.2f us | PW %6.2f us | f %7.3f MHz | DOA %+6.1f deg | %s | purity %.0f %%\n', ...
                k,T.class,T.PRI*1e6,T.pw*1e6,T.freqMed/1e6,T.doa,T.agility,100*T.purity)]; %#ok<AGROW>
        end
    end

% =====================================================================
%                          PLOTTING
% =====================================================================
    function drawAll()
        us  = 1e6;   % s   -> us
        MHz = 1e-6;  % Hz  -> MHz
        P = S.P;
        rr = @(s) real(s);   % complex baseband -> plot the I (real) component

        % --- Tab 1: received ---
        plot(h.ax_rxTime, S.t*us, rr(S.rx), 'Color', [0.2 0.4 0.8]);
        title(h.ax_rxTime, sprintf('Received signal (noisy, SNR %.1f dB, %s)', ...
            P.SNR, P.dom));
        xlabel(h.ax_rxTime, 'Time (\mus)'); ylabel(h.ax_rxTime, 'Amplitude');
        grid(h.ax_rxTime, 'on');

        specPlot(h.ax_rxSpec, S.fFull*MHz, S.magFull, ...
            'Spectrum after filtering', [], S.Fc_est*MHz);

        % --- Tab 2: the filter itself ---
        drawFilter(P);

        % --- Tab 3: band-filtered + envelope ---
        plot(h.ax_filt, S.t*us, rr(S.filt), 'Color', [0.15 0.5 0.35]);
        title(h.ax_filt, sprintf('%s  (%s)', P.filt, P.dom));
        xlabel(h.ax_filt, 'Time (\mus)'); ylabel(h.ax_filt, 'Amplitude');
        grid(h.ax_filt, 'on');

        cla(h.ax_env);
        plot(h.ax_env, S.t*us, S.env, 'Color', [0.6 0.2 0.6]); hold(h.ax_env,'on');
        yline(h.ax_env, S.thr, '--', 'threshold', 'Color', [0.5 0.5 0.5]);
        for k = 1:numel(S.riseAll)
            xr = [S.riseAll(k) S.fallAll(k)]/P.Fs*us;
            area(h.ax_env,xr,[max(S.env) max(S.env)],'FaceAlpha',0.06, ...
                'EdgeColor','none','FaceColor',[0.55 0.55 0.55]);
        end
        for k = 1:numel(S.rise)
            xr = [S.rise(k) S.fall(k)]/P.Fs*us;
            area(h.ax_env,xr,[max(S.env) max(S.env)],'FaceAlpha',0.16, ...
                'EdgeColor','none','FaceColor',[0.9 0.5 0.1]);
        end
        hold(h.ax_env, 'off');
        title(h.ax_env, sprintf('Envelope: %d raw runs, %d selected emitter pulses', S.NdetAll, S.Ndet));
        xlabel(h.ax_env, 'Time (\mus)'); ylabel(h.ax_env, 'Envelope');
        grid(h.ax_env, 'on');

        % --- Tab 4: carrier estimate (zoom to band) ---
        cla(h.ax_freq);
        plot(h.ax_freq, S.fBand*MHz, S.magBand, 'Color', [0.2 0.3 0.7]);
        hold(h.ax_freq, 'on');
        xline(h.ax_freq, S.Fc_est*MHz, '--r', ...
            sprintf('est %.4f MHz', S.Fc_est*MHz), 'LineWidth', 1.2);
        xline(h.ax_freq, (P.Fc+P.fdRx)*MHz, ':k', sprintf('nominal+motion %.4f MHz', (P.Fc+P.fdRx)*MHz));
        hold(h.ax_freq, 'off');
        title(h.ax_freq, sprintf('%s inside search band (sub-bin refined)', P.cest));
        xlabel(h.ax_freq, 'Frequency (MHz)'); ylabel(h.ax_freq, 'Power (dB)');
        grid(h.ax_freq, 'on');

        % --- Tab 5: transmit ---
        GHz = 1e-9;
        plot(h.ax_txTime, S.tx_t*us, real(S.tx), 'Color', [0.75 0.25 0.2]);
        if strcmp(P.txmode, 'Parametric regen')
            title(h.ax_txTime, sprintf(['Parametric regen  |  RF carrier %.4f GHz' ...
                '  (IF %.3f MHz shown)'], S.tx_rf*GHz, P.Ft*MHz));
        else
            title(h.ax_txTime, sprintf(['DRFM replay  |  RF carrier %.4f GHz' ...
                '  (IF %.3f MHz shown, shift %+.3f MHz)'], ...
                S.tx_rf*GHz, P.Ft*MHz, S.tx_shift*MHz));
        end
        xlabel(h.ax_txTime, 'Time (\mus)  [IF / baseband-equivalent]');
        ylabel(h.ax_txTime, 'Amplitude'); grid(h.ax_txTime, 'on');

        % Transmit MIXER output on the RF axis. A single-ended (real-IF) mixer
        % produces BOTH sidebands: multiplying by cos(2*pi*LO*t) copies the IF
        % spectrum to LO +/- f. So the two-sided IF spectrum shifted up by the LO
        % shows the WANTED sideband at f_RF = LO + f_IF and its IMAGE at
        % f_img = LO - f_IF, plus LO leakage at the LO itself. A complex (I/Q)
        % transmit path is one-sided, so its image is naturally absent -- which
        % is exactly why image-reject / IQ upconverters are used.
        cla(h.ax_txSpec);
        [ftw, mtw] = spectrumTwoSided(S.tx, S.P.Fs);
        fRFaxis = (ftw + S.tx_lo) * GHz;              % RF frequency axis (GHz)
        plot(h.ax_txSpec, fRFaxis, mtw, 'Color', [0.25 0.35 0.6]);
        hold(h.ax_txSpec, 'on');
        f_img = (2*S.tx_lo - S.tx_rf);               % image = LO - (RF - LO)
        xline(h.ax_txSpec, S.tx_rf*GHz, '-r', 'wanted', 'LineWidth', 1.2);
        xline(h.ax_txSpec, f_img*GHz, '--', 'image', 'Color', [0.85 0.4 0.1], 'LineWidth', 1.1);
        xline(h.ax_txSpec, S.tx_lo*GHz, ':k', 'LO leak');
        hold(h.ax_txSpec, 'off');
        title(h.ax_txSpec, sprintf('Tx mixer output: wanted %.4f / image %.4f GHz (LO %.4f)', ...
            S.tx_rf*GHz, f_img*GHz, S.tx_lo*GHz));
        xlabel(h.ax_txSpec, 'RF frequency (GHz)'); ylabel(h.ax_txSpec, 'Mag (dB)');
        grid(h.ax_txSpec, 'on'); ylim(h.ax_txSpec, [-80 5]);
        % zoom to a window that shows both sidebands and the LO
        xlim(h.ax_txSpec, [min(f_img, S.tx_rf) - 3*S.P.Fs, ...
                           max(f_img, S.tx_rf) + 3*S.P.Fs]*GHz);

        % --- Tab 6: overview ---
        mini(h.ax_o(1), S.t*us, S.rx,   'Received',     [0.2 0.4 0.8]);
        specPlot(h.ax_o(2), S.fFull*MHz, S.magFull, 'Rx spectrum', ...
            [P.Blo P.Bhi]*MHz, S.Fc_est*MHz);
        mini(h.ax_o(3), S.t*us, S.filt, 'Band-limited', [0.15 0.5 0.35]);
        mini(h.ax_o(4), S.t*us, S.env,  'Envelope',     [0.6 0.2 0.6]);
        plot(h.ax_o(5), S.fBand*MHz, S.magBand, 'Color', [0.2 0.3 0.7]);
        title(h.ax_o(5), 'Carrier est'); grid(h.ax_o(5), 'on');
        xlabel(h.ax_o(5), 'MHz');
        mini(h.ax_o(6), S.tx_t*us, S.tx, 'Transmit',    [0.75 0.25 0.2]);

        % --- Tab 7: geolocation ---
        drawGeo(P);

        % --- Tab 8: time-frequency (spectrograms) ---
        drawSpectrograms(P);

        % --- Tab 9: victim radar range-Doppler payoff ---
        drawVictimRangeDoppler(P);

        % --- Tab 10: receiver LFM pulse compression ---
        drawPulseCompression(P);

        % --- Tab 11: motion, multipath and agility diagnostics ---
        drawEmitterRealism(P);

        % --- Tab 12: complete jammer-versus-ECCM duel ---
        drawEWDuel(P);

        % --- Tab 13: simultaneous-emitter ESM processing ---
        drawMultiEmitterESM(P);
    end

    function drawEmitterRealism(P)
        % Six compact diagnostics turn the realism controls into observable
        % receiver consequences rather than hidden waveform parameters.
        MHz=1e-6; us=1e6;
        ax=h.ax_real;

        cla(ax(1));
        if ~isempty(S.vel_axis)
            plot(ax(1),S.vel_axis,S.doppler_db,'LineWidth',1.1); hold(ax(1),'on');
            xline(ax(1),P.vrad,':k','true');
            xline(ax(1),S.vel_rx_est,'--r','estimated'); hold(ax(1),'off');
            title(ax(1),sprintf('Nonuniform pulse-Doppler: v = %+.1f m/s',S.vel_rx_est));
        else
            title(ax(1),'Need at least two detected pulses for Doppler');
        end
        xlabel(ax(1),'One-way radial velocity (m/s, + closing)'); ylabel(ax(1),'Power (dB)');
        grid(ax(1),'on'); ylim(ax(1),[-60 3]);

        cla(ax(2));
        if ~isempty(S.slow_amp)
            q=1:numel(S.slow_amp);
            ph=unwrap(angle(S.slow_amp))*180/pi;
            m=20*log10(abs(S.mti_amp)/(max(abs(S.slow_amp))+eps)+eps);
            yyaxis(ax(2),'left'); plot(ax(2),q,ph,'-o','MarkerSize',3);
            ylabel(ax(2),'Unwrapped phase (deg)');
            yyaxis(ax(2),'right'); plot(ax(2),q,m,'-x','MarkerSize',4);
            ylabel(ax(2),'2-pulse MTI magnitude (dB)');
            title(ax(2),sprintf('Slow-time phase progression / MTI (%+.1f dB RMS)',S.mti_gain));
        end
        xlabel(ax(2),'Pulse number'); grid(ax(2),'on');

        cla(ax(3));
        ne=numel(S.PRI_seq_est); nt=numel(S.pri_true);
        if nt>0, plot(ax(3),1:nt,S.pri_true*us,'-o','MarkerSize',3); hold(ax(3),'on'); end
        if ne>0, plot(ax(3),1:ne,S.PRI_seq_est*us,'--x','MarkerSize',4); end
        hold(ax(3),'off');
        if nt>0 && ne>0, legend(ax(3),{'true','estimated'},'Location','best');
        elseif nt>0, legend(ax(3),{'true'},'Location','best'); end
        title(ax(3),sprintf('%s: mean %.3f us, std %.3f us',P.primode,S.PRI_est*us,S.PRI_std*us));
        xlabel(ax(3),'Pulse interval'); ylabel(ax(3),'PRI (\mus)'); grid(ax(3),'on');

        cla(ax(4));
        np=numel(S.fcPulse_est);
        if np>0
            plot(ax(4),1:np,S.fcPulse_true*MHz,'-o','MarkerSize',3); hold(ax(4),'on');
            plot(ax(4),1:np,S.fcPulse_est*MHz,'--x','MarkerSize',4); hold(ax(4),'off');
            legend(ax(4),{'true','tracked'},'Location','best');
            title(ax(4),sprintf('%s: tracking RMSE %.1f kHz',P.fhmode,S.carrier_rmse/1e3));
        end
        xlabel(ax(4),'Pulse number'); ylabel(ax(4),'Processing carrier (MHz)'); grid(ax(4),'on');

        cla(ax(5));
        if ~isempty(S.az_perPulse)
            plot(ax(5),1:numel(S.az_perPulse),S.az_perPulse,'-o','MarkerSize',3); hold(ax(5),'on');
            yline(ax(5),P.Az,':k','direct');
            if P.multipath, yline(ax(5),P.mpaz,'--','reflection'); end
            hold(ax(5),'off');
            title(ax(5),sprintf('Per-pulse DF: spread %.2f deg',S.az_spread));
        end
        xlabel(ax(5),'Pulse number'); ylabel(ax(5),'Estimated azimuth (deg)'); grid(ax(5),'on');

        cla(ax(6));
        if ~isempty(S.range_excess)
            plot(ax(6),S.range_excess,S.range_profile_db,'LineWidth',1.1); hold(ax(6),'on');
            xline(ax(6),0,':k','direct');
            if P.multipath, xline(ax(6),P.c*P.mpdelay,'--r','reflection'); end
            hold(ax(6),'off');
            title(ax(6),sprintf('Excess-path profile (multipath loss %.1f dB)',P.mploss));
            xlim(ax(6),[max(min(S.range_excess),-0.25*P.c*P.PW), ...
                min(max(S.range_excess),P.c*(P.mpdelay+1.25*P.PW))]);
        end
        xlabel(ax(6),'Excess one-way path length (m)'); ylabel(ax(6),'Relative magnitude (dB)');
        grid(ax(6),'on'); ylim(ax(6),[-70 3]);
    end

    function drawPulseCompression(P)
        % A direct before/after view of why LFM is useful. The upper axis uses
        % one common time scale so the wide received pulse and narrow matched-
        % filter peak can be compared honestly; the lower axis zooms into the
        % compressed mainlobe and measures its -3 dB width.
        ax1 = h.ax_pcWide; ax2 = h.ax_pcZoom;
        legend(ax1,'off');
        cla(ax1); cla(ax2);

        if ~isfield(S,'pc_active') || ~S.pc_active || ...
                isempty(S.pc_peakTimes) || ~isfinite(S.pc_peakTimes(1))
            text(ax1,0.5,0.55,'Pulse compression is inactive', ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontWeight','bold','FontSize',15);
            text(ax1,0.5,0.43, ...
                'Select "LFM chirp" and choose a non-zero chirp bandwidth.', ...
                'Units','normalized','HorizontalAlignment','center');
            title(ax1,'Before / after matched filtering');
            ax1.XTick=[]; ax1.YTick=[];
            text(ax2,0.5,0.52, ...
                {'A CW rectangular pulse can be integrated, but it has no', ...
                 'LFM time-bandwidth product to create a narrow range spike.'}, ...
                'Units','normalized','HorizontalAlignment','center');
            title(ax2,'Compressed mainlobe');
            ax2.XTick=[]; ax2.YTick=[];
            return
        end

        ax1.XTickMode='auto'; ax1.YTickMode='auto';
        ax2.XTickMode='auto'; ax2.YTickMode='auto';
        us = 1e6;
        t0 = S.pc_peakTimes(1);
        envN = S.env/(max(S.env)+eps);
        pcN = S.pc_mag/(max(S.pc_mag)+eps);
        xWide = [-0.25*P.PW,1.25*P.PW];

        pEnv = plot(ax1,(S.t-t0)*us,envN,'LineWidth',1.2); hold(ax1,'on');
        pPC = plot(ax1,(S.pc_t-t0)*us,pcN,'LineWidth',1.35);
        lPC = xline(ax1,0,'--','compressed TOA');
        legH = [pEnv pPC lPC];
        legT = {'Received envelope','Matched-filter output','Compressed TOA'};
        if isfinite(S.tau_env)
            lEnv = xline(ax1,(S.tau_env-t0)*us,':','envelope TOA');
            legH(end+1) = lEnv; legT{end+1} = 'Envelope TOA';
        end
        hold(ax1,'off');
        xlim(ax1,xWide*us); ylim(ax1,[0 1.08]);
        xlabel(ax1,'Time relative to compressed TOA (\mus)');
        ylabel(ax1,'Normalized magnitude');
        title(ax1,sprintf(['Wide %.2f us LFM pulse -> compressed peak  |  ' ...
            'BT = %.1f, theoretical gain %.1f dB'], ...
            P.PW*us,S.pc_TB,S.pc_gainTheory));
        legend(ax1,legH,legT,'Location','best');
        grid(ax1,'on');

        zoomHalf = max(5/P.bchirp,3*S.pc_width3dB);
        idx = abs(S.pc_t-t0) <= zoomHalf;
        plot(ax2,(S.pc_t(idx)-t0)*us,S.pc_magDB(idx),'LineWidth',1.35);
        hold(ax2,'on');
        yline(ax2,-3,'--','-3 dB');
        xline(ax2,-S.pc_width3dB*us/2,':');
        xline(ax2, S.pc_width3dB*us/2,':');
        hold(ax2,'off');
        xlim(ax2,[-zoomHalf zoomHalf]*us); ylim(ax2,[-45 2]);
        xlabel(ax2,'Time relative to peak (\mus)');
        ylabel(ax2,'Normalized output (dB)');
        title(ax2,sprintf(['Compressed mainlobe: measured -3 dB width %.3f us  |  ' ...
            '1/B = %.3f us'],S.pc_width3dB*us,1/P.bchirp*us));
        grid(ax2,'on');
    end

    function drawVictimRangeDoppler(P)
        % Defended victim map: selected ECCM/MTI output with CFAR detections.
        ax=h.ax_rd; cla(ax); rkm=S.victim_range/1e3; vel=S.victim_velocity;
        imagesc(ax,rkm,vel,S.victim_rdDB.'); set(ax,'YDir','normal');
        try, colormap(ax,'turbo'); catch, colormap(ax,'jet'); end
        caxis(ax,[-50 0]); cb=colorbar(ax); cb.Label.String='Normalized magnitude (dB)';
        hold(ax,'on');
        plot(ax,P.Rng/1e3,0,'wo','MarkerSize',9,'LineWidth',1.5);
        for m=1:size(S.victim_ghostRangePulse,1)
            plot(ax,S.victim_ghostRangePulse(m,:)/1e3,S.victim_ghostVelPulse(m,:), ...
                'w.','MarkerSize',8);
        end
        if ~isempty(S.victim_coverRange)
            plot(ax,S.victim_coverRange/1e3,zeros(size(S.victim_coverRange)),'y+', ...
                'MarkerSize',7,'LineWidth',1.0);
        end
        if ~isempty(S.victim_detRange)
            plot(ax,S.victim_detRange/1e3,S.victim_detVel,'co','MarkerSize',7,'LineWidth',1.2);
        end
        hold(ax,'off');
        xlabel(ax,'Apparent range (km)'); ylabel(ax,'Radial velocity (m/s, + closing)');
        title(ax,sprintf('Defended victim map | ECCM %s | %s | CFAR %s: %d detections', ...
            P.eccm,P.mtiMode,onOff(P.cfar),numel(S.victim_detRange)));
        grid(ax,'on'); if numel(rkm)>1, xlim(ax,[min(rkm) max(rkm)]); end
        if numel(vel)>1, ylim(ax,[min(vel) max(vel)]); end
    end

    function drawEWDuel(P)
        ax=h.ax_duel; rkm=S.victim_range/1e3; vel=S.victim_velocity;
        cla(ax(1)); imagesc(ax(1),rkm,vel,S.victim_rdRawDB.'); set(ax(1),'YDir','normal');
        try, colormap(ax(1),'turbo'); catch, colormap(ax(1),'jet'); end
        caxis(ax(1),[-50 0]); title(ax(1),'Raw victim map: clutter + ECM');
        xlabel(ax(1),'Range (km)'); ylabel(ax(1),'Velocity (m/s)'); grid(ax(1),'on');

        cla(ax(2)); imagesc(ax(2),rkm,vel,S.victim_rdDB.'); set(ax(2),'YDir','normal');
        try, colormap(ax(2),'turbo'); catch, colormap(ax(2),'jet'); end
        caxis(ax(2),[-50 0]); hold(ax(2),'on');
        if ~isempty(S.victim_detRange)
            plot(ax(2),S.victim_detRange/1e3,S.victim_detVel,'co','MarkerSize',6,'LineWidth',1.1);
        end
        plot(ax(2),P.Rng/1e3,0,'wo','MarkerSize',8,'LineWidth',1.3); hold(ax(2),'off');
        if P.cfar
            title(ax(2),sprintf('Defended + CFAR: skin %s, ghosts %d/%d', ...
                hitStr(S.victim_skinDetected),S.victim_ghostDetected,numel(S.victim_ghostRange)));
        else
            title(ax(2),'Defended map (CFAR disabled)');
        end
        xlabel(ax(2),'Range (km)'); ylabel(ax(2),'Velocity (m/s)'); grid(ax(2),'on');

        cla(ax(3)); p=1:numel(S.victim_trackPeak);
        plot(ax(3),p,S.victim_trackPeak/1e3,'--','LineWidth',1.0); hold(ax(3),'on');
        plot(ax(3),p,S.victim_trackLead/1e3,'LineWidth',1.3);
        plot(ax(3),p,S.victim_trackUsed/1e3,'k.','MarkerSize',12);
        yline(ax(3),P.Rng/1e3,':','skin'); hold(ax(3),'off');
        legend(ax(3),{'Peak tracker','Leading edge','Selected','True skin'},'Location','best');
        xlabel(ax(3),'Pulse number'); ylabel(ax(3),'Tracked range (km)');
        title(ax(3),sprintf('RGPO tracker duel: %s selected',trackerName(P))); grid(ax(3),'on');

        cla(ax(4));
        yyaxis(ax(4),'left');
        if ~isempty(S.victim_ghostRangePulse)
            plot(ax(4),p,S.victim_ghostRangePulse(1,:)/1e3,'-o','MarkerSize',3);
        else
            plot(ax(4),p,nan(size(p)));
        end
        ylabel(ax(4),'Ghost range (km)');
        yyaxis(ax(4),'right'); plot(ax(4),p,S.victim_vgpoVel,'-s','MarkerSize',3);
        ylabel(ax(4),'VGPO velocity (m/s)'); xlabel(ax(4),'Pulse number');
        title(ax(4),sprintf('PRI/frequency ECCM: %d ambiguous track(s) rejected',S.victim_priReject));
        grid(ax(4),'on');
    end

    function s=hitStr(tf)
        if isfield(S,'P') && ~S.P.cfar
            s='CFAR OFF';
        elseif tf
            s='DETECTED';
        else
            s='not detected';
        end
    end

    function drawSpectrograms(P)
        % Short-time Fourier transform of the received capture (top) and the
        % retransmitted signal (bottom). With an LFM chirp you can SEE the
        % diagonal sweep; DRFM replays it (shifted in frequency, delayed in
        % time, possibly several false targets), while parametric regen shows a
        % flat CW line -- a direct picture of what DRFM preserves.
        MHz = 1e-6; us = 1e6;
        if P.iq
            Wc=0.5*(P.Bhi-P.Blo); rxBand=[-Wc Wc]*MHz;
        else
            rxBand=[P.Blo P.Bhi]*MHz;
        end
        drawOneSpectrogram(h.ax_sgRx,S.filt,S.t,P,'Received capture',rxBand);
        drawOneSpectrogram(h.ax_sgTx, S.tx, S.tx_t, P, ...
            sprintf('Retransmit IF/baseband-equiv (%s, RF %.4f GHz)', ...
                P.txmode, S.tx_rf*1e-9), []);
    end

    function drawGeo(P)
        % LEFT: a radar-scope (PPI-style) half-plane showing the TRUE source
        % position and the ESTIMATED one in (range, azimuth) space. Azimuth is
        % measured from boresight (straight up = 0 deg), +/-90 deg to the sides.
        us = 1e6;
        ax = h.ax_ppi; cla(ax); hold(ax, 'on');

        Rt = P.Rng/1e3;                       % true range (km)
        Re = S.Rng_est/1e3;                   % est  range (km)
        azt = P.Az; aze = S.Az_est;

        % DRFM FALSE TARGETS -- read off the ACTUAL transmitted waveform. The
        % transmit chain ran the pulse detector on S.tx and stored the detected
        % apparent ranges in S.tx_ghostR (km after the /1e3 below). Because these
        % come from real detection, overlapping/interfering false targets MERGE
        % here exactly as they do on the wire, instead of showing ideal delays.
        ghostR = [];
        isDRFM = ~strcmp(P.txmode, 'Parametric regen');
        if isDRFM && isfinite(aze) && isfield(S,'tx_ghostR') && ~isempty(S.tx_ghostR)
            gr = S.tx_ghostR / 1e3;                 % m -> km (apparent range)
            % show them only when there is actual deception (delay or >1 target)
            if P.ntgt > 1 || P.drdelay > 0
                ghostR = gr(:).';
            end
        end

        cand = [Rt, Re, ghostR]; cand = cand(isfinite(cand) & cand > 0);
        if isempty(cand), Rmax = 1; else, Rmax = 1.25 * max(cand); end

        % range rings (semicircles) + labels
        thArc = linspace(-90, 90, 121);
        for rr = Rmax*(0.25:0.25:1)
            plot(ax, rr*sind(thArc), rr*cosd(thArc), '-', ...
                'Color', [0.85 0.85 0.85]);
            text(ax, 0, rr, sprintf('%.1f km', rr), 'Color', [0.5 0.5 0.5], ...
                'FontSize', 8, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
        end
        % azimuth spokes every 30 deg
        for a = -90:30:90
            plot(ax, [0 Rmax*sind(a)], [0 Rmax*cosd(a)], '-', ...
                'Color', [0.9 0.9 0.9]);
            text(ax, 1.05*Rmax*sind(a), 1.05*Rmax*cosd(a), sprintf('%d', a), ...
                'Color', [0.5 0.5 0.5], 'FontSize', 8, ...
                'HorizontalAlignment', 'center');
        end
        plot(ax, 0, 0, '^', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.2 0.2], ...
            'MarkerEdgeColor', [0.2 0.2 0.2]);   % the receive array (origin)

        % true source (circle) and estimate (x), plus a bearing line to the est
        if isfinite(Rt)
            plot(ax, Rt*sind(azt), Rt*cosd(azt), 'o', 'MarkerSize', 11, ...
                'LineWidth', 1.5, 'Color', [0.15 0.45 0.85]);
            text(ax, Rt*sind(azt), Rt*cosd(azt), '  true', ...
                'Color', [0.15 0.45 0.85], 'FontSize', 9);
        end
        if isfinite(Re) && isfinite(aze)
            plot(ax, [0 Re*sind(aze)], [0 Re*cosd(aze)], '-', ...
                'Color', [0.85 0.30 0.20]);
            plot(ax, Re*sind(aze), Re*cosd(aze), 'x', 'MarkerSize', 13, ...
                'LineWidth', 2, 'Color', [0.85 0.30 0.20]);
            text(ax, Re*sind(aze), Re*cosd(aze), '  est', ...
                'Color', [0.85 0.30 0.20], 'FontSize', 9);
        end

        % DRFM false targets (ghosts) as magenta diamonds down the bearing.
        for g = 1:numel(ghostR)
            gx = ghostR(g)*sind(aze); gy = ghostR(g)*cosd(aze);
            plot(ax, gx, gy, 'd', 'MarkerSize', 9, 'LineWidth', 1.3, ...
                'Color', [0.75 0.1 0.75], 'MarkerFaceColor', [0.95 0.75 0.95]);
            text(ax, gx, gy, sprintf('  ghost %d', g), ...
                'Color', [0.6 0.1 0.6], 'FontSize', 8);
        end

        hold(ax, 'off');
        axis(ax, 'equal');
        xlim(ax, [-1.15*Rmax 1.15*Rmax]); ylim(ax, [-0.1*Rmax 1.2*Rmax]);
        ax.XTick = []; ax.YTick = [];
        if isempty(ghostR)
            title(ax, sprintf('Source location  (R err %+.2f km, Az err %+.2f deg)', ...
                (S.Rng_est-P.Rng)/1e3, S.Az_est-P.Az));
        else
            title(ax, sprintf('Source + %d DRFM false target(s)  (Az err %+.2f deg)', ...
                numel(ghostR), S.Az_est-P.Az));
        end

        % RIGHT: the two array channels overlaid on the first detected pulse,
        % so the inter-element phase offset that drives the azimuth estimate is
        % actually visible.
        ax2 = h.ax_chan; cla(ax2); hold(ax2, 'on');
        if S.Ndet >= 1
            i1 = S.rise(1); i2 = min(numel(S.t), S.fall(1));
            % show a short window around the first pulse's leading portion
            span = i1:min(i2, i1 + max(20, round(3*P.Fs/max(S.Fc_est,eps))));
            plot(ax2, S.t(span)*us, S.filt(span),  'Color', [0.15 0.45 0.85]);
            plot(ax2, S.t(span)*us, S.filt2(span), 'Color', [0.85 0.30 0.20]);
            legend(ax2, {'ch1', 'ch2'}, 'Location', 'northeast');
            title(ax2, sprintf('Array channels @ first pulse  (\\Delta\\phi = %.1f deg)', ...
                rad2deg(S.dphi_est)));
        else
            title(ax2, 'No pulse detected -> no azimuth');
        end
        hold(ax2, 'off');
        xlabel(ax2, 'Time (\mus)'); ylabel(ax2, 'Amplitude'); grid(ax2, 'on');
    end

    function drawFilter(P)
        % Show the EFFECTIVE magnitude response of whichever filter is selected
        % (effective = what the signal actually experiences, including the
        % |H|^2 squaring from forward-backward IIR filtering) plus, where one
        % exists, the impulse response.
        MHz = 1e-6;
        f = linspace(0, P.Fs/2, 2048);
        [Hmag, himp] = filterResponse(P, f);
        HdB = 20*log10(Hmag + 1e-12);

        cla(h.ax_hmag);
        plot(h.ax_hmag, f*MHz, HdB, 'Color', [0.2 0.2 0.6], 'LineWidth', 1.2);
        hold(h.ax_hmag, 'on');
        xline(h.ax_hmag, P.Blo*MHz, '--', 'Color', [0.5 0.5 0.5]);
        xline(h.ax_hmag, P.Bhi*MHz, '--', 'Color', [0.5 0.5 0.5]);
        yline(h.ax_hmag, -3, ':', '-3 dB', 'Color', [0.7 0.3 0.3]);
        xline(h.ax_hmag, P.Fc*MHz, ':k', 'carrier');
        hold(h.ax_hmag, 'off');
        title(h.ax_hmag, sprintf('%s  |  effective magnitude response', P.filt));
        xlabel(h.ax_hmag, 'Frequency (MHz)'); ylabel(h.ax_hmag, 'Mag (dB)');
        ylim(h.ax_hmag, [-100 10]); grid(h.ax_hmag, 'on');

        cla(h.ax_himp);
        if isempty(himp)
            % FFT-mask filters have no finite impulse response worth drawing --
            % say so rather than showing a misleading empty axis.
            title(h.ax_himp, 'FFT-mask filter: applied in the frequency domain');
            text(h.ax_himp, 0.5, 0.5, ...
                {'This filter multiplies the FFT of the entire record.', ...
                 'It is non-causal and has no realisable FIR/IIR kernel.', ...
                 'Switch to FIR or Butterworth to see an impulse response.'}, ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
            h.ax_himp.XTick = []; h.ax_himp.YTick = [];
        else
            n = 0:numel(himp)-1;
            stem(h.ax_himp, n, himp, 'filled', 'MarkerSize', 2, ...
                'Color', [0.15 0.5 0.35]);
            title(h.ax_himp, sprintf('Impulse response (%d samples)', numel(himp)));
            xlabel(h.ax_himp, 'Sample'); ylabel(h.ax_himp, 'h[n]');
            grid(h.ax_himp, 'on');
        end
    end

    function updateResults()
        MHz=1e-6; us=1e6; P=S.P;
        ftruth=P.Fc+P.fdRx;
        ferr=(S.fc_center_est-ftruth)*MHz;
        if isnan(S.PRI_est), prf=NaN; else, prf=1/S.PRI_est; end
        if isnan(S.PW_est)||isnan(S.PRI_est), duty=NaN;
        else, duty=100*S.PW_est/S.PRI_est; end
        pwErr=(S.PW_est-P.PW)*us; priErr=(S.PRI_est-P.PRI)*us;
        rngErr=(S.Rng_est-P.Rng)/1e3; azErr=S.Az_est-P.Az;
        mpTxt='OFF'; if P.multipath, mpTxt=sprintf('ON, %.1f dB loss, %.2f us / %.1f deg', ...
                P.mploss,P.mpdelay*us,P.mpaz); end

        txt=sprintf([ ...
            'FILTER : %s\nCARRIER EST : %s\nWAVEFORM : %s\n\n' ...
            'MEASURED (Rx)\n-------------------------\n' ...
            'Global peak : %8.4f MHz\n' ...
            'Pulse centre: %8.4f MHz\n' ...
            'Expected ctr: %8.4f MHz (motion included)\n' ...
            'Centre error: %+8.4f MHz\n' ...
            'Carrier track RMSE: %8.2f kHz\n' ...
            'Pulse runs  : %8d raw / %d selected (true %d)\n' ...
            'Pulse width : %8.3f us (err %+.3f)\n' ...
            'PRI mean    : %8.3f us (err %+.3f)\n' ...
            'PRI std     : %8.3f us\n' ...
            'PRF         : %8.3f kHz\nDuty cycle  : %8.2f %%\n\n' ...
            'EMITTER REALISM\n-------------------------\n' ...
            'Radial vel  : true %+.2f / est %+.2f m/s\n' ...
            'Rx Doppler  : true %+.3f / est %+.3f kHz (one-way)\n' ...
            'PRI agility : %s, depth %.1f %%\n' ...
            'Freq agility: %s, step %.3f MHz\n' ...
            'Multipath   : %s\n' ...
            'DF pulse spread: %8.3f deg\n' ...
            'MTI RMS ratio  : %8.2f dB\n\n' ...
            'GEOLOCATION\n-------------------------\n' ...
            'Range est   : %8.3f km (err %+.3f)\n' ...
            '   (true     %8.3f km)\nRange method: %s\n' ...
            'Envelope rng: %8.3f km (err %+.3f)\n' ...
            'Azimuth est : %8.2f deg (err %+.2f)\n' ...
            '   (true     %8.2f deg)\n' ...
            'Baseline    : %8.2f m  (%.2f lambda)\n' ...
            'DF cal      : %s\n   uncal az : %8.2f deg (err %+.2f)\n\n' ...
            'PULSE COMPRESSION  (main Rx)\n-------------------------\n' ...
            'Status      : %s\nPulses scored: %8d\nTime-bandwd : %8.2f\n' ...
            'Theory gain : %8.2f dB\nInput SNR est : %8.2f dB\n' ...
            'Output SNR est: %8.2f dB\nMeasured gain : %8.2f dB\n' ...
            '-3 dB width : %8.3f us\nTOA res 1way: %8.2f m\nRadar res   : %8.2f m\n\n' ...
            'TRANSMIT  (IF -> RF upconversion)\n-------------------------\n' ...
            'Mode        : %s\nTx IF       : %8.3f MHz\nDRFM shift  : %8.3f MHz\n' ...
            'Upconv LO   : %10.5f GHz\nTx RF carr  : %10.5f GHz\n' ...
            'RGPO delay  : %8.3f us\nDoppler     : %8.3f kHz\nFalse tgts  : %8d\n\n' ...
            'VICTIM RADAR  (matched filter + Doppler FFT)\n-------------------------\n' ...
            'Range res   : %8.2f m\nVelocity res: %8.2f m/s\n' ...
            'Unambig vel : +/-%7.2f m/s\nGhost returns: %8d\nGhost vel   : %8.2f m/s'], ...
            P.filt,P.cest,waveformStr(P),S.Fc_est*MHz,S.fc_center_est*MHz, ...
            ftruth*MHz,ferr,S.carrier_rmse/1e3,S.NdetAll,S.Ndet,P.Np, ...
            S.PW_est*us,pwErr,S.PRI_est*us,priErr,S.PRI_std*us,prf/1e3,duty, ...
            P.vrad,S.vel_rx_est,P.fdRx/1e3,S.fd_rx_est/1e3,P.primode,100*P.pridev, ...
            P.fhmode,P.hopstep*MHz,mpTxt,S.az_spread,S.mti_gain, ...
            S.Rng_est/1e3,rngErr,P.Rng/1e3,S.rangeMethod,S.Rng_env/1e3, ...
            (S.Rng_env-P.Rng)/1e3,S.Az_est,azErr,P.Az,P.base, ...
            P.base/(P.c/max(P.Flo+S.Fc_ctr,eps)),calStatusStr(P),S.Az_raw,S.Az_raw-P.Az, ...
            pulseCompressionStatus(),S.pc_Ncoh,S.pc_TB,S.pc_gainTheory,S.pc_snrIn, ...
            S.pc_snrOut,S.pc_gainMeasured,S.pc_width3dB*us,S.pc_toaRes,S.pc_radarRes, ...
            P.txmode,P.Ft*MHz,S.tx_shift*MHz,S.tx_lo*1e-9,S.tx_rf*1e-9, ...
            P.drdelay*us,P.dop/1e3,round(P.ntgt),S.victim_rangeRes,S.victim_velRes, ...
            S.victim_velUnamb,numel(S.victim_ghostRange),victimGhostVelocity());
        duel=sprintf([ ...
            '\n\nEW DUEL  (victim ECCM / ECM)\n-------------------------\n' ...
            'ECCM mode    : %s\nVictim PRI dev: %7.2f %%\nVictim hop   : %8.3f MHz\n' ...
            'Ground clutter: %s (%+.1f dB, %.1f m/s spread)\nClutter cancel: %s\n' ...
            'CFAR         : %s, Pfa %.1e\nCFAR cells   : %8d\nClutter suppression: %7.2f dB\n' ...
            'Skin detected: %s\nGhosts detected: %d / %d\nPRI-ambig tracks: %d\nPRI cells rejected: %d\n' ...
            'VGPO slope   : %8.3f kHz/pulse\nCover pulses : %s (%d, %.1f dB loss)\n' ...
            'Noise jammer : %s, J/S %+.1f dB, spot BW %.3f MHz\n' ...
            'Tracker      : %s\nMedian tracked range: %.3f km'], ...
            P.eccm,100*P.vpriDev,P.vhopstep/1e6,onOff(P.clutter),P.clutterDB,P.clutterV, ...
            P.mtiMode,onOff(P.cfar),P.pfa,numel(S.victim_detRange),S.victim_clutterSuppression,hitStr(S.victim_skinDetected), ...
            S.victim_ghostDetected,numel(S.victim_ghostRange),S.victim_priReject,S.victim_priCellReject, ...
            P.vgstep/1e3,onOff(P.cover),P.coverN,P.coverLoss,P.jammode,P.jsdb,P.spotbw/1e6, ...
            trackerName(P),median(S.victim_trackUsed)/1e3);
        txt=[txt,duel,esmResultsText()];
        if ~isempty(S.notes)
            txt=[txt,sprintf('\n\nNOTES\n-------------------------\n')];
            for k=1:numel(S.notes), txt=[txt,sprintf('- %s\n',S.notes{k})]; end %#ok<AGROW>
        end
        h.results.Text=txt;
    end

    function s = pulseCompressionStatus()
        if isfield(S,'pc_active') && S.pc_active
            s = 'ACTIVE (LFM matched filter)';
        elseif strcmp(S.P.modtype,'LFM chirp')
            s = 'INACTIVE (no valid compressed peak)';
        else
            s = 'OFF (select LFM chirp)';
        end
    end

    function v = victimGhostVelocity()
        if isfield(S,'victim_ghostVel') && ~isempty(S.victim_ghostVel)
            v = S.victim_ghostVel(1);
        else
            v = NaN;
        end
    end

    function s = calStatusStr(P)
        % Human-readable calibration status for the results panel.
        if ~P.usecal
            s = sprintf('OFF  (ch phase offset %.0f deg)', P.phoff*180/pi);
        elseif cal.valid
            s = sprintf('ON   (fine %+.1f / short %+.1f deg removed)', ...
                cal.dphiF*180/pi, cal.dphiS*180/pi);
        else
            s = 'ON but NOT calibrated -> press CALIBRATE';
        end
    end

% =====================================================================
%              FILTER BANK  (all hand-rolled, no toolbox)
% =====================================================================
end
