% Finding Segment Boundaries Function
function boundaries = findSegmentBoundaries(x, y, window)
    % find_segment_boundaries() defines the bounds of each segment boundary
    % given an osx-CAAX-GFP signal along a bony ray by using the second
    % derivative of the signal.
    % Inputs: y (column vector), the signal.
    %         window (double), the size of the window over which to search
    %         for a bound starting from the middle of the segment boundary.
    % Outputs: boundaries (column vector), contains all the segment
    % boundaries in the form [left bound of segment boundary 1, right bound
    % of segment boundary 1, left bound of segment boundary 2, ... right
    % bound of segment boundary n].
    smooth_y = smoothdata(y, 'gaussian'); % Applies a Gaussian filter to the signal.
    inv_smooth_y = -smooth_y;
    d2_smooth_y = smoothdata(del2(smooth_y), 'gaussian', 100); % Applies a Gaussian filter to the second derivative.
    inv_d2_smooth_y = -d2_smooth_y; % Inversion of the second derivative to find minima.

    min_prom = 5;
    min_width = 200;
    [~, idxs] = findpeaks(inv_smooth_y, x, 'MinPeakProminence', min_prom, 'MinPeakDistance', min_width); % Returns vector of local maximums w/ prominence > alpha.
    boundaries = zeros(height(idxs)*2, 1);
    %beta = 0.001;
    for j = 1:height(idxs) % Iterates over all the maxima.
        idx = idxs(j, :);
        if idx-window > 1
            [~, x1, ~, ~] = findpeaks(inv_d2_smooth_y(idx-window:idx, :), 'NPeaks', 1); % Finds minimum before the maximum.
        else
            [~, x1, ~, ~] = findpeaks(inv_d2_smooth_y(1:idx, :), 'NPeaks', 1); % Finds minimum before the maximum.
        end
        if idx+window > height(inv_d2_smooth_y)
            [~, x2, ~, ~] = findpeaks(inv_d2_smooth_y(idx:height(inv_d2_smooth_y), :), 'NPeaks', 1); % Finds minimum after the maximum.
        else
            [~, x2, ~, ~] = findpeaks(inv_d2_smooth_y(idx:idx+window, :), 'NPeaks', 1); % Finds minimum after the maximum.
        end
        x1 = idx - window + x1; % Normalizes to the correct index (not the windown index).
        x2 = idx + x2; % Normalizes to the correct index (not the windown index).
        boundaries((j*2)-1, 1) = x1;
        boundaries(j*2, 1) = x2;
    end
end