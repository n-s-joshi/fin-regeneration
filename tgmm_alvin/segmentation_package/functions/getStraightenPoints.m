function [IM, points] = getStraightenPoints(IM)
% This function takes in an image and gets the points required for
% straightening by straighten(). Written by NSJ on 03142025, adapted from 
% Jiacheng Wang.
    a = 200; % How much you want to buffer the IM matrix by. This is to prevent the hitting the edge of the image during straightening.
    horizontal_buffer = zeros(height(IM), a);
    vertical_buffer = zeros(a, numel(IM(1, :)) + 2*a);
    IM = [horizontal_buffer, IM, horizontal_buffer];
    IM = [vertical_buffer; IM; vertical_buffer];
    f = figure(1);clf;
    f.WindowState = 'maximized';
    imagesc(IM);axis image off;
    hold on;
    [x,y,b]=ginput(1);
    plot(x,y,'o-');hold off
    while b~=3 % Left click on all points except the last one. Right click the last one. Then hit 'enter/return'.
        [x_,y_,b]=ginput(1);
        if b==8
            x(end)=[];y(end)=[];
        else
            x=[x,x_];y=[y,y_];
        end
        f.Children(1).Children(1).XData=x;
        f.Children(1).Children(1).YData=y;
    end
    points = [x;y];
end