function [N,fstop,Hdb,f,b,w,info] = find_ntaps(fc,fstop,fs,Adb,win)
% Usage: [N,fstop,Hdb,f,b,w,info] = find_ntaps(fc,fstop,fs,Adb,win)
%
% Iteratively search for a windowed sinc FIR filter satisfying
% provided specifications with the fewest number of coefficients.
%
%   fc...........cutoff frequency spec in Hz, the frequency at
%                which the normalized magnitude response is -6 dB
%   fstop........stopband frequency in Hz, fstop is the frequency at
%                which the normalized magnitude response is down
%                Adb db or more
%   fs...........sampling frequency Hz
%   Adb..........stopband attenuation in dB
%   win..........char array containing built-in window name
%                built-in windows: 'hann', 'hamming','kaiser',
%                'blackman','parzen','nuttall','albrecht'
%                alternatively, a char array containing a function 
%                handle such as '@tukeywin'
%   N............output: number of coefficients
%   fstop........output: pass-through of fstop input
%   Hdb..........output: 1-d array magnitude response
%   f............output: 1-d array of frequencies for which the
%                magnitude reponse is evaluated
%   b............output: 1-d array of filter coefficients
%   w............output: 1-d array of window coefficients
%   info.........output: struct of useful info
%

    if nargin ~= 5
        fprintf(2,'Wrong number of input arguments.\n\n');
        help find_ntaps
        return
    end

    if fstop <= fc
        fprintf(2,'fstop must be greater than fc.\n');
        return
    end

    if 2*fstop >= fs
        fprintf(2,'fstop must be in the range (fc,fs/2)\n');
        return
    end

    if exist('freqz') == 2
        fresp = @freqz;
    else
        fresp = @pm_freqz;
    end

    fc0 = fc / fs;
    fstop0 = fstop / fs;
    tbw = fstop0 - fc0;
    tbw2 = tbw / (2*fc0);

    N = ceil(1/tbw * (Adb - 8) / 14); % fharris eq. 3.11
    if rem(N,2) == 0
        N = N + 1;
    end
    N = max(N,11); % ensure min and odd

    nfreq = 1000;
    f0 = fstop0 + (-nfreq:nfreq)/nfreq * tbw;
    f0(nfreq + 1) = fstop0;
    f0 = f0(f0>0 & f0<0.5);

    iters = 0;
    tol = 1 / 100;
    %incr = ceil(N/2);
    incr = N;
    if rem(incr,2) == 1
        incr = incr + 1;
    end
    best_err = inf;
    best_N   = N;
    best_ix = nan;
    early_term = 0;
    err = 1;
    sgn = 1;
    while 1
        iters = iters + 1;

        [b,w,N] = wsinc(fc0,1,Adb,win,N);
        Hdb = 20*log10(abs(fresp(b,1,f0,1)));
        ix = find(Hdb > -Adb, 1, 'last');
        if isempty(ix)
            N = 2 * N + iters;
            if rem(N,2) == 0
                N = N + 1;
            end
        elseif ix == numel(Hdb)
            N = ceil(1.1 * N);
            if rem(N,2) == 0
                N = N + 1;
            end
        else
            tbw_cand = (f0(ix) - fc0)/(2*fc0);
            %err = tbw_cand / tbw2 - 1;
            err = f0(ix) - fstop0;
            if abs(err) < abs(best_err)
                if isnan(best_ix)
                    sgn = sign(err);
                end
                best_err = err;
                best_N = N;
                best_ix = ix;
            end
            %if inrange(err, -tol, 0)
            if inrange(err, -1e-9, 0)
                break
            elseif err > 0
                if sgn == -1
                    incr = max(ceil(incr/2),2);
                end
                if rem(incr,2) == 1
                    incr = incr + 1;
                end
                N = N + incr;
            else % err < 0
                if sgn == 1
                    incr = max(ceil(incr/2),2);
                end
                while incr >= N
                    incr = max(ceil(incr/2),2);
                end
                if rem(incr,2) == 1
                    incr = incr - 1;
                end
                incr = max(incr,2);
                N = N - incr;
                N = max(N,11);
            end
            if err > 0
                sgn = 1;
            else
                sgn = -1;
            end
        end
        if N > 100000
            % Too many coefficients. Design did not converge.
            break
        end
        if iters >= 50
            if best_err > 0
                N = best_N + 1;
            else
                N = best_N;
            end
            err = best_err;
            ix = best_ix;
            [b,w,nused,wsincinfo] = wsinc(fc0,1,Adb,win,N);
            break
        end

        nvec(iters) = N;
        errvec(iters) = err;

        % Terminate early if stuck in a limit cycle
        if iters > 4
            nhist = nvec(iters-3:end) - nvec(iters-3);
            if all(nhist == [0 2 0 2]) ...
            || all(nhist == [0 0 0 0]) ...
            || all(nhist == [0 -1 0 -1])
                early_term = 1;
                break
            end
        end
           
    end % while

    [b,w,N,wsincinfo] = wsinc(fc0,1,Adb,win,N);
    Hdb = 20*log10(abs(fresp(b,1,f0,1)));

    if Hdb(nfreq+1) + Adb > 1
        fprintf(2,'Failed to converge (win=%s, err=%.4f, N=%d Hdb(fstop)=%.1f).\n', ...
            win,err,N,Hdb(nfreq+1));
    end

    % collect some useful info
    info = struct();
    info.win = win;
    info.winbeta = wsincinfo.beta;
    info.nvec = nvec;
    info.errvec = errvec;

    % Final frequency response
    f = (0:2^18)/2^18 * fs/2;
    Hdb = 20*log10(abs(fresp(b,1,f,fs)));

end % function
