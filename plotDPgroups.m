function [] = plotDPgroups(f2cand, f2Idx, fRatio, fs)

dpCoeffs = [
   -1  1;
    2 -1;    % 2*f1 - f2
    1  0;    % f1
    0  1;    % f2
   -1  2;
    2  0;
    1  1;
    0  2;
    3  0;
    2  1;
    1  2;
    0  3
];

% ---------- Build f1 and f2 pairs ----------
pairs  = [f2cand./fRatio, f2cand];
nPairs = size(pairs,1);

nyq      = fs/2;
% binWidth = fs/blockLength;
% Nbins    = blockLength/2;

nDPs     = size(dpCoeffs,1);
freqsAll = nan(nPairs,nDPs);

for i = 1:nPairs
    freqs = dpCoeffs(:,1)*pairs(i,1) + dpCoeffs(:,2)*pairs(i,2);
    freqs(freqs<=0 | freqs>=nyq) = NaN;
    freqsAll(i,:) = freqs.';
end

k = size(f2Idx, 1);
groups = cell(k,1);

for s = 1:k
    groupVec = f2Idx(s,:)';
    groupVec = groupVec(groupVec~=0);
    groups{s} = groupVec;
end


% ---------- Simple visualization (highlights main freqs f1,f2,DP=2f1-f2) ----------
figure('Name','Scheduled groups (frequencies)');
hold on;

% Slot colors
nSlots = max(1,numel(groups));
colors = lines(nSlots);

% Marker shapes for roles
markerMain = {'s','d','^'};         % markers for f1,f2,DP
roleNames = {'f1','f2','2f1-f2'};

% visual params
alphaFactor = 0.45;        % how faded the "other" frequencies are (0..1)
mainMarkerSize = 9;
otherMarkerSize = 6;
mainLineWidth = 1.2;

% Build an ordering that sorts pairs by their assigned slot.
% Pairs in the same slot keep their original index order.
slotOfPair = zeros(size(pairs,1),1);
for s = 1:numel(groups)
    slotOfPair(groups{s}) = s;
end
[~, sortIdx] = sortrows([slotOfPair, (1:size(pairs,1))']); % primary by slot, secondary by original index
sortedPairs = sortIdx(:);

% We'll map original pair index -> y-position (1..nPairs) according to sorting by slot
yPos = nan(size(pairs,1),1);
for y = 1:numel(sortedPairs)
    yPos(sortedPairs(y)) = y;
end

for s = 1:numel(groups)
    idxs = groups{s}(:);
    baseColor = colors(mod(s-1,size(colors,1))+1,:);
    fadedColor = (1-alphaFactor)*baseColor + alphaFactor*[1 1 1]; % fade toward white

    for k = 1:numel(idxs)
        i = idxs(k);
        fRow = freqsAll(i, :);
        valid = ~isnan(fRow);
        fRow = fRow(valid);

        % y position determined by slot-sorted ordering
        y = yPos(i);

        % identify main frequencies for this pair
        f1 = pairs(i,1);
        f2 = pairs(i,2);
        dp = 2*f1 - f2;

        % Plot "other" frequencies (all valid except the three mains) as faded markers
        isMain = abs(fRow - f1) < 1e-9 | abs(fRow - f2) < 1e-9 | abs(fRow - dp) < 1e-9;
        others = fRow(~isMain);
        if ~isempty(others)
            plot(others, y*ones(size(others)), 'o', ...
                'Color', fadedColor, 'MarkerFaceColor', fadedColor, 'MarkerSize', otherMarkerSize);
        end

        % Plot main frequencies with distinct shapes but colored by slot
        plot(f1, y, markerMain{1}, 'MarkerSize', mainMarkerSize, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor', baseColor, 'LineWidth', mainLineWidth);
        plot(f2, y, markerMain{2}, 'MarkerSize', mainMarkerSize, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor', baseColor, 'LineWidth', mainLineWidth);
        plot(dp, y,  markerMain{3}, 'MarkerSize', mainMarkerSize, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor', baseColor, 'LineWidth', mainLineWidth);

        % Connect all frequencies for this pair with a thin line in slot color
        if numel(fRow) > 1
            fRowSorted = sort(fRow);
            plot(fRowSorted, y*ones(size(fRowSorted)), '-', 'Color', baseColor, 'LineWidth', 0.8);
        end
    end
end

set(gca,'xscale','log');
xlabel('Frequency (Hz)');
ylabel('Frequency pair (sorted by slot)'); % indicate y-axis meaning
yticks(1:size(pairs,1));
yticklabels(arrayfun(@num2str, sortedPairs, 'UniformOutput', false)); % show original pair indices in sorted order
grid on; box on;

% Legend:
% - slot color patches
% - role markers
legendHandles = gobjects(0);
legendLabels = {};

% slot color entries (show up to first 10 slots to avoid huge legends)
maxLegendSlots = min(nSlots, 10);
for s = 1:maxLegendSlots
    h = plot(nan, nan, 's', 'MarkerFaceColor', colors(s,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    legendHandles(end+1) = h; %#ok<AGROW>
    legendLabels{end+1} = sprintf('Slot %d', s); %#ok<AGROW>
end

% role markers
for r = 1:numel(markerMain)
    h = plot(nan, nan, markerMain{r}, 'MarkerSize', mainMarkerSize, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', [0 0 0]); % black face for role legend
    legendHandles(end+1) = h; %#ok<AGROW>
    legendLabels{end+1} = roleNames{r}; %#ok<AGROW>
end

legend(legendHandles, legendLabels, 'Location', 'bestoutside');