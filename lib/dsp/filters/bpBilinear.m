    function [b, a] = bpBilinear(pLP, n, f1, f2, Fs)
        % Shared analog-lowpass-prototype -> digital-bandpass machinery:
        %   1. prewarp the band edges (bilinear warps frequency)
        %   2. lowpass -> bandpass transform: each prototype pole p becomes two,
        %      solving  s^2 - p*BW*s + W0^2 = 0
        %   3. bilinear transform poles/zeros into the z-plane
        % Zeros: n at s=0 (from the BP transform) -> z=+1, and n at s=inf -> z=-1.
        W1 = 2*Fs*tan(pi*f1/Fs);
        W2 = 2*Fs*tan(pi*f2/Fs);
        BW   = W2 - W1;
        W0sq = W1 * W2;

        pb = zeros(1, 2*n);
        for ii = 1:n
            half = pLP(ii)*BW/2;
            disc = sqrt(half^2 - W0sq);
            pb(2*ii-1) = half + disc;
            pb(2*ii  ) = half - disc;
        end
        zb = zeros(1, n);

        fs2 = 2*Fs;
        pz = (1 + pb/fs2) ./ (1 - pb/fs2);
        zz = (1 + zb/fs2) ./ (1 - zb/fs2);       % s = 0   -> z = +1
        zz = [zz, -ones(1, n)];                  % s = inf -> z = -1

        b = real(poly(zz));
        a = real(poly(pz));
    end
