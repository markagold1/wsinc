function p = plotres(res,figno)

    if nargin == 2
        figure(figno)
    else
        figure
    end
    
    f = res{1}.f;
    Hdb = nan(numel(f),numel(res));
    leg = {};
    for kk = 1:numel(res)
       Hdb(:,kk) = res{kk}.Hdb(:);
       leg{end+1} = res{kk}.window;
    end
    p = plot(f,Hdb,'LineWidth',1.5);
    hleg = legend(leg);
    if ~isoctave()
        xtickformat('usd');
        xtickformat('%.0f');
    end
    grid on
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
    title('Magnitude Responses');
    set(gca,'FontSize',14);
    set(hleg,'FontSize',10);

end % function
