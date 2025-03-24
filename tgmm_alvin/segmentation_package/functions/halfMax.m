function index = halfMax(y)
    [value, max_index] = max(y);
    half_max = value*0.6;
    residuals = (abs(y(1:max_index, 1) - half_max));
    [~, index] = min(residuals);
end