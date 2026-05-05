% This script should be replaced by frequency distribution algorithm
f2aim = [500:500:1500 2000:1000:10000].';
% f2Idx = [2 3 4 6 7 8 11 12; 1 5 9 10 zeros(1,4)];
f2Idx = [2 4 6 10; 3 5 9 0; 7 12 0 0];
save("dpFrequencies.mat", "f2aim", "f2Idx");