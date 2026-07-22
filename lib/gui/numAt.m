    function fld = numAt(grid, row, labCol, fldCol, txt, val, fmt)
        % Labelled numeric field placed at explicit label/field columns.
        if nargin < 7, fmt = '%.6g'; end
        l = uilabel(grid, 'Text', txt);
        l.Layout.Row = row; l.Layout.Column = labCol;
        fld = uieditfield(grid, 'numeric', 'Value', val, ...
            'ValueDisplayFormat', fmt);
        fld.Layout.Row = row; fld.Layout.Column = fldCol;
    end
