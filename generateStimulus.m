function [stimulusOut] = generateStimulus(nFreq, freq, t)
    stimulus = [0 0];
    for i = nFreq
        stimulus = stimulus + [cos(2*pi*freq(i,1)*t) cos(2*pi*freq(i,2)*t)];
    end
    
    stimulusOut = stimulus/i;
end