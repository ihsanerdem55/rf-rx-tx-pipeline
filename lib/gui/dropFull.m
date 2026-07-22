    function dd = dropFull(grid, row, txt, items, val)
        % Dropdown given the full width (label col 1, menu spans cols 2-4) so
        % long filter/estimator names are never truncated.
        l = uilabel(grid, 'Text', txt);
        l.Layout.Row = row; l.Layout.Column = 1;
        dd = uidropdown(grid, 'Items', items, 'Value', val);
        dd.Layout.Row = row; dd.Layout.Column = [2 4];
    end
