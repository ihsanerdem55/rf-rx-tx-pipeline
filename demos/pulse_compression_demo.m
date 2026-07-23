function pulse_compression_demo(doPlots)
% PULSE_COMPRESSION_DEMO  Headless demonstration of pulse compression.
%
%   Two teaching experiments, printed as quantitative tables (and, if a
%   display is available, drawn as compressed-pulse profiles):
%
%     PART A  Receive weighting on an LFM matched filter -- the range
%             sidelobe / mainlobe-width / SNR trade-off. The window is applied
%             to the reference only (a mismatched filter), exactly as the GUI
%             app does, so the achievable floor is limited to about -30 dB.
%
%     PART B  Range sidelobes by waveform: CW (no compression) vs LFM
%             (~ -14 dB sidelobes) vs a Barker-13 phase code (~ -20 dB).
%
%   Runs in base MATLAB with NO GUI and NO Signal Processing Toolbox, so it
%   can be used to validate the DSP core headlessly. It reuses the app's own
%   taper, lib/dsp/pulse/pcWindow.m.
%
%   USAGE (from the repository root):
%       >> addpath demos
%       >> pulse_compression_demo          % tables (+ plots if a display exists)
%       >> pulse_compression_demo(false)   % tables only

    if nargin < 1, doPlots = true; end
    thisDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(thisDir, '..', 'lib')));

    Fs = 40e6; PW = 10e-6; B = 3e6;         % 40 MHz Fs, 10 us pulse, 3 MHz sweep
    L  = round(PW*Fs);
    tl = (0:L-1)/Fs - PW/2;
    lfm = exp(1i*pi*(B/PW)*tl.^2);          % baseband LFM reference

    % ---- PART A: receive weighting on the LFM matched filter --------------
    windows = {'None','Taylor','Hann','Hamming','Blackman'};
    fprintf('\nPART A  Receive weighting on an LFM matched filter (TB = %.0f)\n', PW*B);
    fprintf('        signal = plain LFM, filter = windowed LFM (mismatched)\n');
    fprintf('  %-9s %10s %12s %12s\n', 'window', 'PSL (dB)', '-3dB (us)', 'broadening');
    w0 = NaN;
    for k = 1:numel(windows)
        refW = lfm .* pcWindow(L, windows{k});
        mag  = abs(conv(lfm, conj(fliplr(refW))));
        [psl, wSamp] = localMetrics(mag);
        if k == 1, w0 = wSamp; end
        fprintf('  %-9s %10.1f %12.3f %12s\n', windows{k}, psl, ...
            wSamp/Fs*1e6, sprintf('x%.2f', wSamp/max(w0,1)));
    end

    % ---- PART B: waveform comparison (matched filter, no weighting) -------
    fprintf('\nPART B  Compressed range sidelobes by waveform (matched filter)\n');
    fprintf('  %-12s %8s %10s %12s\n', 'waveform', 'TB', 'PSL (dB)', '-3dB (us)');
    waves = buildWaveforms(Fs, PW, B, L);
    for k = 1:numel(waves)
        mag = abs(conv(waves(k).x, conj(fliplr(waves(k).x))));
        [psl, wSamp] = localMetrics(mag);
        fprintf('  %-12s %8.0f %10s %12.3f\n', waves(k).name, waves(k).TB, ...
            pslStr(psl), wSamp/Fs*1e6);
    end
    fprintf('\n');

    if doPlots
        try
            plotCompression(Fs, L, lfm, windows, waves);
        catch err
            fprintf('(plots skipped -- no display? %s)\n', err.message);
        end
    end
end

% =====================================================================
function [pslDB, width] = localMetrics(mag)
% Compressed-pulse metrics: -3 dB mainlobe width (samples) and peak sidelobe
% level (dB). PSL is measured BEYOND the first nulls, not the -3 dB edges --
% the skirt down to the first null is part of the mainlobe.
    mag = mag(:).';
    [pk, ip] = max(mag);
    thr = pk/sqrt(2);
    il = ip; while il > 1        && mag(il-1) >= thr, il = il-1; end
    ir = ip; while ir < numel(mag) && mag(ir+1) >= thr, ir = ir+1; end
    width = ir - il + 1;
    nl = ip; while nl > 1        && mag(nl-1) <= mag(nl), nl = nl-1; end
    nr = ip; while nr < numel(mag) && mag(nr+1) <= mag(nr), nr = nr+1; end
    side = mag([1:nl-1, nr+1:end]);
    if isempty(side)
        pslDB = NaN;                        % e.g. CW: triangular autocorrelation
    else
        pslDB = 20*log10(max(side)/max(pk,eps));
    end
end

% =====================================================================
function waves = buildWaveforms(Fs, PW, B, L)
    tl  = (0:L-1)/Fs - PW/2;
    cw  = ones(1, L);                                   % rectangular CW, TB = 1
    lfm = exp(1i*pi*(B/PW)*tl.^2);                      % LFM chirp
    Pn   = struct('modtype', 'NLFM chirp',      'bchirp', B, 'PW', PW);
    Pb   = struct('modtype', 'Barker-13 code',  'bchirp', 0, 'PW', PW);
    nlfm = exp(1i*chirpPhase(tl, Pn));                  % NLFM (shared lib design)
    code = exp(1i*chirpPhase(tl, Pb));                  % Barker-13 (shared lib code)
    waves = struct( ...
        'name', {'CW (rect)', 'LFM chirp', 'NLFM chirp', 'Barker-13'}, ...
        'x',    {complex(cw), lfm, nlfm, code}, ...
        'TB',   {1, PW*B, PW*B, 13});
end

% =====================================================================
function s = pslStr(psl)
    if isfinite(psl), s = sprintf('%.1f', psl); else, s = 'n/a'; end
end

% =====================================================================
function plotCompression(Fs, L, lfm, windows, waves)
    us = 1e6;
    figure('Name', 'Pulse compression demo', 'Color', 'w');

    subplot(1, 2, 1); hold on;
    for k = 1:numel(windows)
        refW = lfm .* pcWindow(L, windows{k});
        mag  = abs(conv(lfm, conj(fliplr(refW))));
        mag  = mag / max(mag);
        t = ((1:numel(mag)) - (numel(mag)+1)/2) / Fs * us;
        plot(t, 20*log10(mag+eps), 'LineWidth', 1.1, 'DisplayName', windows{k});
    end
    hold off; grid on; xlim([-4 4]); ylim([-70 2]);
    xlabel('Time (\mus)'); ylabel('Normalized (dB)');
    title('PART A: LFM receive weighting'); legend('Location', 'best');

    subplot(1, 2, 2); hold on;
    for k = 1:numel(waves)
        mag = abs(conv(waves(k).x, conj(fliplr(waves(k).x))));
        mag = mag / max(mag);
        t = ((1:numel(mag)) - (numel(mag)+1)/2) / Fs * us;
        plot(t, 20*log10(mag+eps), 'LineWidth', 1.1, 'DisplayName', waves(k).name);
    end
    hold off; grid on; xlim([-12 12]); ylim([-40 2]);
    xlabel('Time (\mus)'); ylabel('Normalized (dB)');
    title('PART B: waveform sidelobes'); legend('Location', 'best');
end
