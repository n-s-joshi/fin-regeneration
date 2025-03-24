function [IM2, w] = optimizeStraighten(IM, points, w)
% This function find a good width to input into the straighten function if
% you do not know a decent width a priori. Written by NSJ, 03172025.
    optimal = false;
    i = 1;
    while (optimal == false) && (i < 100)
        % disp(['Iterations: ' num2str(i)]);
        IM2 = straighten(IM, points, w);
        if (mean(IM2(:, 1)) < 1) && (mean(IM2(:, end)) < 1)
            optimal = true;
        end
        w = w + 10;
        i = i + 1;
    end
end