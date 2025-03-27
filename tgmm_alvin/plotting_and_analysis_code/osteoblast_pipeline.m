%% Prepare
clear;clc;close all;
disp('------------Preparing------------');

% Path to directory containing the script
addpath(genpath('/Volumes/NSJ_Data_I/scripts/tgmm_alvin/segmentation_package/functions/'));
disp('Done');

% Init plot settings
set(groot,'defaultLineLineWidth',2.0)
num_times = 6;
newcolors = brewermap(num_times, '+Reds');

x_conversion = 0.606; % Conversion from pixels to um.

%% Load Paths
disp('------------Loading Paths------------');
% Init dictionary to hold paths
paths =[];

% Base path
paths.diskName = '/Volumes/NSJ_Data_I';
paths.expName = 'caudal_fin/11092024_osx-caax-gfp_reamp';
paths.objFolder = [paths.diskName filesep paths.expName filesep 'objects'];
paths.boundariesFolder = [paths.diskName filesep paths.expName filesep 'boundaries'];
paths.alignedFolder = [paths.diskName filesep paths.expName filesep 'aligned'];

% Paths to  new directories
paths.csvFolder = [paths.diskName filesep paths.expName filesep 'csv'];
paths.plotsFolder = [paths.diskName filesep paths.expName filesep 'plots'];

% Make directories for those which do not already exist
mkdir(paths.plotsFolder);
mkdir(paths.boundariesFolder);
mkdir(paths.alignedFolder);

if isfile([paths.objFolder filesep 'analysis_mat.mat'])
    load([paths.objFolder filesep 'analysis_mat.mat']);
end

disp('Done');

%% Create structure to hold data
disp('------------Initializing Data------------');

% Naturally sort .mat files 
st_dir = dir([paths.objFolder filesep 'fish*_ray*_*hpa.mat']);
[~,sort_idx] = sort_nat({st_dir.name});
st_dir = st_dir(sort_idx);

% Structure which will hold all data
analysis_mat = [];

for i = 1:numel(st_dir)
    % print
    disp('---------------')
    disp(st_dir(i).name)

    analysis_mat_here = [];

    nameHere = st_dir(i).name;
    frt = sscanf(st_dir(i).name,'fish%d_ray%d_%dhpa'); % fish ray time
    analysis_mat_here.fish = frt(1);
    analysis_mat_here.ray = frt(2);
    analysis_mat_here.hpa = frt(3);
    analysis_mat_here.dpa = frt(3)./24;

    load([paths.objFolder filesep 'fish' num2str(frt(1)) '_ray' num2str(frt(2)) '_' num2str(frt(3)) 'hpa' '.mat']); %load myScale

    analysis_mat_here.amp_plane = myScale.manualAmpPlane;
    analysis_mat_here.straightened_x_length = myScale.straightenedXLength;
    analysis_mat_here.x_pixels = (1:height(myScale.linescan))';
    analysis_mat_here.x_microns = analysis_mat_here.x_pixels * x_conversion;
    analysis_mat_here.raw_profile = myScale.linescan;

    analysis_mat = [analysis_mat analysis_mat_here];
end
save([paths.objFolder filesep 'analysis_mat'], 'analysis_mat');
disp('Done');

%% Align profiles
disp('------------Aligning Profiles by Amputation Plane------------')
j = 1;
f = figure;
colororder(newcolors);
for i = 1:width(analysis_mat)
    name = ['fish' num2str(analysis_mat(i).fish) '_ray' num2str(analysis_mat(i).ray) '_' num2str(analysis_mat(i).hpa) 'hpa'];
    disp(name);
    
    if analysis_mat(i).hpa == 0
        reference_amp_x = analysis_mat(i).straightened_x_length - analysis_mat(i).amp_plane(1, 1);
    end
    
    amp_x = analysis_mat(i).straightened_x_length - analysis_mat(i).amp_plane(1, 1);
    shift = amp_x - reference_amp_x;
    shift = round(shift);
    analysis_mat(i).shift = shift;
    analysis_mat(i).x_pixels_shifted = analysis_mat(i).x_pixels - shift;
    analysis_mat(i).x_microns_shifted = analysis_mat(i).x_pixels_shifted * x_conversion;
    
    if j < num_times
        plot(analysis_mat(i).x_microns_shifted, analysis_mat(i).raw_profile); hold on;
        j = j + 1;
    elseif j == num_times
        plot(analysis_mat(i).x_microns_shifted, analysis_mat(i).raw_profile); hold off;
        title(['fish' num2str(analysis_mat(i).fish) '_ray' num2str(analysis_mat(i).ray) '_aligned'], 'Interpreter', 'none');
        saveas(f, [paths.alignedFolder filesep 'fish' num2str(analysis_mat(i).fish) '_ray' num2str(analysis_mat(i).ray) '_aligned.png']);
        close(f);
        j = 1;
        f = figure;
        colororder(newcolors);
    end

end
save([paths.objFolder filesep 'analysis_mat'], 'analysis_mat');

disp('Done');
    
%% Find and Exclude Segment Boundaries
disp('------------Identifying and Excluding Segment Boundaries------------');
% j = 1;
for k = 1:width(analysis_mat)
    name = ['fish' num2str(analysis_mat(k).fish) '_ray' num2str(analysis_mat(k).ray) '_' num2str(analysis_mat(k).hpa) 'hpa'];
    disp(name);
    window = 100;
    if analysis_mat(k).hpa == 0 % Only finds boundaries on profiles from 0hpa and propagates to all the subsequent timepoints.
        boundaries = findSegmentBoundaries(analysis_mat(k).x_pixels, analysis_mat(k).raw_profile, window);
    end
    
    boundaries_shifted = boundaries + analysis_mat(k).shift;
    new_y = analysis_mat(k).raw_profile(1:boundaries_shifted(1));
    for i = 1:height(boundaries_shifted)/2
        x1 = boundaries_shifted(i*2, :);
        if i == height(boundaries_shifted)/2
            x2 = height(analysis_mat(k).x_pixels);
        else
            x2 = boundaries_shifted((i*2)+1, :);
        end
        new_y = [new_y; NaN(x1-height(new_y)-1, 1); analysis_mat(k).raw_profile(x1:x2)];
    end
    smooth_y = feval(fit_spline(analysis_mat(k).x_pixels, new_y), analysis_mat(k).x_pixels);

    analysis_mat(k).smooth_excluded_y = smooth_y;

    f = figure;
    figure(f);
    plot(analysis_mat(k).x_pixels, analysis_mat(k).raw_profile);
    hold on
    plot(analysis_mat(k).x_pixels, analysis_mat(k).smooth_excluded_y);
    for bounds = boundaries_shifted(:, 1) %*x_conversion
        xline(bounds);
        hold on
    end
    hold off
    title(name, 'Interpreter', 'none');
    xlim([1, 1200]);
    ylim([1, 255]);
    saveas(f, [paths.boundariesFolder filesep name '_boundaries.png']);
    close;
    % j = j+1;
end
save([paths.objFolder filesep 'analysis_mat'], 'analysis_mat');
disp('Done');

%% Plot Smoothened Ray Profiles
disp('------------Plotting Smoothened Profiles------------');
video = VideoWriter([paths.plotsFolder filesep 'fish' num2str(analysis_mat(1).fish) '_ray' num2str(analysis_mat(1).ray) '_video.mp4'], 'MPEG-4');
video.FrameRate = 1;
video.Quality = 100;
open(video);
f = figure;
colororder(newcolors);
for i = 1:width(analysis_mat)
    plot(analysis_mat(i).x_microns_shifted, analysis_mat(i).smooth_excluded_y);
    xlim([-400, 750]);
    ylim([0, 255]);
    title(['fish' num2str(analysis_mat(i).fish) '_ray' num2str(analysis_mat(i).ray) '_' num2str(analysis_mat(i).dpa) 'dpa'], 'Interpreter', 'none');
    xlabel('Distance from amputation plane (microns)');
    ylabel('Osx:CAAX-GFP intensity (A.U.)');
    legend('0hpa', '12hpa', '24hpa', '48hpa', '60hpa', '72hpa');
    frame = getframe(f); %get frame
    writeVideo(video, frame);
    if j < num_times
        hold on;
        j = j + 1;
    elseif j == num_times
        hold off;
        close(f);
        close(video);
        j = 1;
        video = VideoWriter([paths.plotsFolder filesep 'fish' num2str(analysis_mat(i+1).fish) '_ray' num2str(analysis_mat(i+1).ray) '_video.mp4'], 'MPEG-4');
        video.FrameRate = 1;
        video.Quality = 100;
        open(video);
        f = figure;
        colororder(newcolors);
    end
end
disp('Done');

%% Half-max Quantification and Plotting
newcolors = ["red"; "red"; "blue"; "blue"];

x = zeros(width(timepoints), 1);
for i=1:width(timepoints)
    t = timepoints{1, i};
    t = str2num(t);
    x(i, 1) = t;
end

j = 1;
for fish = fish_nums
    f = figure;
    colororder(newcolors);
    g = figure;
    colororder(newcolors);
    h = figure;
    colororder(newcolors);
    for ray = ray_nums
        for j = 1:num_bifurcates(j)
            raw_y = zeros(width(timepoints), 1);
            smoothed_y = zeros(width(timepoints), 1);
            smoothed_excluded_y = zeros(width(timepoints), 1);
            for k = 1:width(timepoints)
                y_name = strjoin(['fish' fish '_' 'ray' ray '_' num2str(j) '_' timepoints{1, k} 'hpa' '_y'], '');
                raw_y(k, 1) = halfMax(raw_table{:, y_name})*x_conversion;
                smoothed_y(k, 1) = halfMax(smoothed_table{:, y_name})*x_conversion;
                smoothed_excluded_y(k, 1) = halfMax(smoothed_excluded_table{:, y_name})*x_conversion;
            end
            figure(f);
            plot(x, raw_y);
            hold on
            figure(g);
            plot(x, smoothed_y);
            hold on
            figure(h);
            plot(x, smoothed_excluded_y);
            hold on
        end
        j = j + 1;
    end
    figure(f);
    legend('Ray 2', 'Ray 2', 'Ray 6', 'Ray 6');
    xticks(x);
    xlabel('Time (hpa)');
    ylabel('Inflection Point (microns away from amputation plane)');
    title(strjoin(['fish' fish ' ' 'Raw'], ''));
    figure(g);
    legend('Ray 2', 'Ray 2', 'Ray 6', 'Ray 6');
    xticks(x);
    xlabel('Time (hpa)');
    ylabel('Inflection Point (microns away from amputation plane)');
    title(strjoin(['fish' fish ' ' 'Smoothened'], ''));
    figure(h);
    legend('Ray 2', 'Ray 2', 'Ray 6', 'Ray 6');
    xticks(x);
    xlabel('Time (hpa)');
    ylabel('Inflection Point (microns away from amputation plane)');
    title(strjoin(['fish' fish ' ' 'Smoothened and Excluded'], ''));
end
disp('Done');

%% Test
x = zeros(width(timepoints), 1);
for i=1:width(timepoints)
    t = timepoints{1, i};
    t = str2num(t);
    x(i, 1) = t;
end

j = 1;
for fish = fish_nums
    for ray = ray_nums
        for j = 1:num_bifurcates(j)
            for k = 1:width(timepoints)
                x_name = strjoin(['fish' fish '_' 'ray' ray '_' num2str(j) '_' timepoints{1, k} 'hpa' '_x'], '');
                y_name = strjoin(['fish' fish '_' 'ray' ray '_' num2str(j) '_' timepoints{1, k} 'hpa' '_y'], '');
                
                figure();
                plot(raw_table{:, x_name}, raw_table{:, y_name}, 'Color', 'red');
                hold on
                xline(halfMax(raw_table{:, y_name}), 'red');
                hold on
                plot(smoothed_table{:, x_name}, smoothed_table{:, y_name}, 'Color', 'blue');
                hold on
                xline(halfMax(smoothed_table{:, y_name}), 'blue');
                hold on
                plot(smoothed_excluded_table{:, x_name}, smoothed_excluded_table{:, y_name}, 'Color', 'green');
                hold on
                xline(halfMax(smoothed_excluded_table{:, y_name}), 'green');
            end
        end
        j = j + 1;
    end
end
disp('Done');