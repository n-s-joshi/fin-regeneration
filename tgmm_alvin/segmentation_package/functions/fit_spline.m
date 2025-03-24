function output = fit_spline(x, y)
    % spline() fits a smoothing spline to a given signal.
    % Inputs: x (column vector), x range of a signal.
    %         y (column vector), the signal.
    % Outputs: output (cfit object), the function of the spline fit.
    [xData, yData] = prepareCurveData(x, y);
    
    % Set up fittype and options.
    ft = fittype('smoothingspline');
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 0.0000005; % Inversely proportional to how much weight should be given to the second derivative ("curviness").
    
    % Fit model to data.
    [fitresult, gof] = fit(xData, yData, ft, opts);
    output = fitresult;
end
