function [] = simpleFFTPlot(signal, fs, xlims)
    if nargin < 3; clear x1 x2; end

    n = length(signal);
    Y = abs(fft(signal, n));
    X = fs*(1:(n/2))/n;
    
    if nargin == 3
        semilogx(X, 20*log10(abs(Y(1:floor(n/2)))))
        xlim(xlims);
    else
        semilogx(X, 20*log10(abs(Y(1:n/2))))
    end
end