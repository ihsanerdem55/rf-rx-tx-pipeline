    function sec4(grid, row, txt)
        % Bold section header spanning all four control columns.
        l = uilabel(grid, 'Text', txt, 'FontWeight', 'bold', ...
            'FontColor', [0.15 0.30 0.55]);
        l.Layout.Row = row; l.Layout.Column = [1 4];
    end
