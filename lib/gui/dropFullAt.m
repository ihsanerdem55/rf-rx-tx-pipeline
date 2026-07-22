    function dd = dropFullAt(grid, row, labCol, fldCol, txt, items, val)
        % Compact dropdown occupying one label/field pair.
        l = uilabel(grid, 'Text', txt);
        l.Layout.Row = row; l.Layout.Column = labCol;
        dd = uidropdown(grid, 'Items', items, 'Value', val);
        dd.Layout.Row = row; dd.Layout.Column = fldCol;
    end
