function p = plotres(res,varargin)

    if nargin == 1
        figure;
    end

    xscale = 'linear';
    for kk = 1:numel(varargin)
        thisarg = varargin{kk};
        if isnumeric(thisarg)
            figure(thisarg);
        elseif ischar(thisarg) && ...
          (strcmpi(thisarg,'linear') || strcmpi(thisarg,'log'))
            xscale = lower(thisarg);
        end
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
    if strcmpi(xscale,'log')
        set(gca,'XScale','log');
    end

end % function
