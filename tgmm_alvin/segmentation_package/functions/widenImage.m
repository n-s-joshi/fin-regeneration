function [IM2] = widenImage(IM)
    a = 200; % How much you want to buffer the IM matrix by. This is to prevent the hitting the edge of the image during straightening.
    horizontal_buffer = zeros(height(IM), a);
    vertical_buffer = zeros(a, numel(IM(1, :)) + 2*a);
    IM2 = [horizontal_buffer, IM, horizontal_buffer];
    IM2 = [vertical_buffer; IM2; vertical_buffer];
end