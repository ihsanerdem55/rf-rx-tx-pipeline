    function s = sincLocal(x)
        % sinc(x) = sin(pi x)/(pi x), with the removable singularity at x = 0.
        % (MATLAB's sinc() lives in the Signal Processing Toolbox.)
        s = ones(size(x));
        nz = (x ~= 0);
        s(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
    end
