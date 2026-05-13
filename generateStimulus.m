function [stimulusOut] = generateStimulus(nFreq, freq, t, levelPerTone)
    stimulus = [0 0];
    if nargin <= 3
        for i = nFreq
            stimulus = stimulus + [cos(2*pi*freq(i,1)*t) cos(2*pi*freq(i,2)*t)];
        end
    else
        levels = levelPerTone.*sqrt(2);
        for i = nFreq
            stimulus = stimulus + [cos(2*pi*freq(i,1)*t).*levels(i, 1) cos(2*pi*freq(i,2)*t).*levels(i, 2)];
        end
    end
    
    stimulusOut = stimulus/i;
end