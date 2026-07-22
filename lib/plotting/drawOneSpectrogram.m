    function drawOneSpectrogram(ax, x, tvec, P, ttl, bandMHz)
        us = 1e6; MHz = 1e-6;
        [Tc, Fc_ax, Sdb] = stftLocal(x, P.Fs);
        cla(ax);
        imagesc(ax, Tc*us, Fc_ax*MHz, Sdb);
        set(ax, 'YDir', 'normal');
        try, colormap(ax, 'turbo'); catch, colormap(ax, 'jet'); end
        if ~isempty(bandMHz)
            hold(ax, 'on');
            yline(ax, bandMHz(1), '--w'); yline(ax, bandMHz(2), '--w');
            hold(ax, 'off');
        end
        if isreal(x), ylim(ax,[0 P.Fs/2*MHz]);
        else, ylim(ax,[-P.Fs/2 P.Fs/2]*MHz); end
        title(ax, ttl);
        xlabel(ax, 'Time (\mus)'); ylabel(ax, 'Frequency (MHz)');
    end
