function g = scanGain(t, P)
% SCANGAIN  Antenna-scan amplitude modulation imprinted on the emitted pulse
%   train as the emitter's beam sweeps relative to the receiver. t is a scalar
%   or vector of emission times (s); g is the per-time gain (0..~1.6). Base
%   MATLAB only.
%
%     'Circular'  rotating search beam: one Gaussian mainlobe flash per
%                 revolution (period 1/scanrate), flash width set by beamwidth.
%     'Conical'   nutating tracker: sinusoidal amplitude modulation at the scan
%                 (nutation) rate -- the modulation a conical-scan tracker uses
%                 for angle sensing, and that ECM exploits.
%     'Off'       g = 1 (no modulation; default -> unchanged pulse train).
    switch P.scanmode
        case 'Circular'
            frac = mod(P.scanrate*t, 1);            % position within a revolution
            off  = 360*min(frac, 1-frac);           % degrees from Rx boresight
            sig  = max(P.scanbw*0.42466, eps);      % Gaussian sigma from -3 dB beamwidth
            g = exp(-0.5*(off./sig).^2);
        case 'Conical'
            g = max(0, 1 + 0.6*cos(2*pi*P.scanrate*t));
        otherwise
            g = ones(size(t));
    end
    g = reshape(g, size(t));
end
