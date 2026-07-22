    function specPlot(ax, f, magdB, ttl, bandMHz, markMHz)
        cla(ax);
        plot(ax, f, magdB, 'Color', [0.25 0.35 0.6]); hold(ax, 'on');
        if ~isempty(bandMHz)
            xline(ax, bandMHz(1), '--', 'Color', [0.5 0.5 0.5]);
            xline(ax, bandMHz(2), '--', 'Color', [0.5 0.5 0.5]);
        end
        if ~isempty(markMHz) && ~isnan(markMHz)
            xline(ax, markMHz, '-r', 'LineWidth', 1.1);
        end
        hold(ax, 'off');
        title(ax, ttl);
        xlabel(ax, 'Frequency (MHz)'); ylabel(ax, 'Mag (dB)');
        grid(ax, 'on'); ylim(ax, [-80 5]);
    end
