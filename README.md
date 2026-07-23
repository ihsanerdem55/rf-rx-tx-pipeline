# rf-rx-tx-pipeline

Interactive RF receive / retransmit simulation — an educational DSP + EW demo
built entirely on base MATLAB (no Signal Processing Toolbox). It walks a pulsed
RF signal through a complete pipeline (generate → receive/DF → retransmit →
victim radar → multi-emitter ESM) and draws graphics at every stage.

Requires **MATLAB R2021a or later** (`uifigure`, scrollable `uigridlayout`,
`uitabgroup`).

## Running

```matlab
>> rf_rxtx_pipeline
```

The entry point adds `lib/` (and its subfolders) to the path automatically, so
you only need `rf_rxtx_pipeline.m` on your path or to be in this
folder when you launch it.

## Demos

Standalone scripts under `demos/` run in **base MATLAB with no GUI and no
display**, so the DSP core can be exercised and validated headlessly:

```matlab
>> addpath demos
>> pulse_compression_demo          % prints tables (+ plots if a display exists)
>> pulse_compression_demo(false)   % tables only
```

`pulse_compression_demo` shows two teaching experiments — receive weighting on
an LFM matched filter (the range-sidelobe / mainlobe-width trade-off, reusing
the app's `pcWindow` taper) and a waveform comparison (CW vs LFM vs NLFM vs a
Barker-13 phase code). The **NLFM chirp** is also selectable in the app's *Mod
type* control; its stationary-phase design lowers range sidelobes (~ -24 dB)
without the SNR loss of a receive window.

## Project layout

The app used to be a single ~3800-line file. It is now split into a small app
core plus a library of single-purpose functions, one function per file, grouped
by concern:

```
rf_rxtx_pipeline.m                   App core: GUI construction, the pipeline
                                    driver, the four processing "chains", and
                                    all drawing/results code. These functions
                                    share the GUI state (model S, view h,
                                    calibration cal) as nested functions, so
                                    they stay together in one file.

lib/
  gui/        uifigure/uigridlayout builder helpers (labels, numeric fields,
              dropdowns, paired-axes) used to lay out the control panel.
  util/       small shared helpers: input validation, struct merge, waveform /
              on-off / tracker-name formatting.
  signal/     ground-truth signal generation: pulse-train synthesis, emitter
              PRI/hop schedule, azimuth resolution.
  dsp/
    filters/  from-scratch band-pass filters and their machinery — FFT
              brick-wall / cosine-taper / Gaussian masks, FIR windowed-sinc /
              Kaiser, Butterworth & Chebyshev-I IIR (analog prototype →
              bilinear), plus freqz/impz/filtfilt-lite and sinc.
    spectral/ carrier estimation and spectra: periodogram & Welch PSD,
              one/two-sided spectrum, analytic signal.
    pulse/    envelope pulse detection, main-train selection, and LFM
              pulse-compression (match filtering).
  tx/         transmit-side jamming aids: cover pulses, noise jammer.
  victim/     coherent pulse-Doppler victim radar: PRI build, clutter, CA-CFAR,
              range trackers, slow-time DFT, PRI-consistency weighting.
  esm/        multi-emitter ESM: emitter profiles, scene simulation, PDW
              extraction, deinterleaving (k-means), track stats/classification,
              MUSIC DoA, STFT.
  plotting/   reusable plot helpers (spectrogram tiles, spec/mini plots).
```

## How the split works

MATLAB nested functions share the parent workspace and must live in the same
file, so the functions that read or write the shared GUI state (`S`, `h`,
`cal`) remain nested inside `rf_rxtx_pipeline.m`. Every function
that is **pure** — it takes its inputs as arguments and returns its outputs,
touching none of the shared state — was moved out to its own file under `lib/`.
The extraction was verified to be behaviour-preserving: each moved function's
code is byte-for-byte identical to the original, and the call graph is closed
(no library function calls back into the app core).

> Note: the `lib/` functions deliberately reimplement several Signal Processing
> Toolbox routines under plain names (`spectrum`, `analytic`, `filtfiltLite`,
> …). They are added to the path only while the app is on the path; this is the
> same "base-MATLAB-only" design the original single file relied on.
