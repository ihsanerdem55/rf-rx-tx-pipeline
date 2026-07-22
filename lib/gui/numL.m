    function fld = numL(grid, row, txt, val, fmt)   % left pair  (cols 1,2)
        if nargin < 5, fmt = '%.6g'; end
        fld = numAt(grid, row, 1, 2, txt, val, fmt);
    end
