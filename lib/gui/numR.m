    function fld = numR(grid, row, txt, val, fmt)   % right pair (cols 3,4)
        if nargin < 5, fmt = '%.6g'; end
        fld = numAt(grid, row, 3, 4, txt, val, fmt);
    end
