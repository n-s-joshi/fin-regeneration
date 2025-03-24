%% Prepare
clear;clc;close all;
disp('------------Preparing------------');

% Path to directory containing the script
addpath(genpath('/Volumes/NSJ_Data_I/scripts'));
disp('Done');

% Init plot settings
set(groot,'defaultLineLineWidth',2.0)

x_conversion = 0.606; % Conversion from pixels to um.

%% Load Paths
disp('------------Loading Paths------------');
% Init dictionary to hold paths
paths =[];

% Base path
paths.diskName = '/Volumes/NSJ_Data_I';
paths.expName = 'caudal_fin/11092024_osx-caax-gfp_reamp';
paths.objFolder = [paths.diskName filesep paths.expName filesep 'objects'];

% Paths to  new directories
paths.csvFolder = [paths.diskName filesep paths.expName filesep 'csv'];
paths.plotsFolder = [paths.diskName filesep paths.expName filesep 'plots'];

% Make directories for those which do not already exist
mkdir(paths.plotsFolder);
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
    analysis_mat_here.x_pixels = (1:height(myScale.linescan))';
    analysis_mat_here.x_microns = analysis_mat_here.x_pixels * x_conversion;
    analysis_mat_here.raw_profile = myScale.linescan;

    analysis_mat = [analysis_mat analysis_mat_here];
end
% 
% 
% fish_nums = {'1', '2', '3'};
% ray_nums = {'2', '6'};
% % timepoints = {'0', '4', '8', '12', '16', '20', '24'}; 
% timepoints = {'0', '12', '24', '48', '60', '72'};
% %num_bifurcates = [2, 3, 2, 2, 2, 1]; 
% num_bifurcates = [2, 2, 2, 2, 2, 2]; % How many 'subrays' are being quantified in each image?
% 
% var_names_raw_data = {};
% var_names = {};
% i = 1;
% for fish = fish_nums
%     for ray = ray_nums
%         for j = 1:num_bifurcates(i)
%             for time = timepoints
%                 name = ['fish' fish '_ray' ray '_' num2str(j) '_' time 'hpa'];
%                 new_column_x = [name '_x'];
%                 new_column_y = [name '_y'];
%                 var_names{end+1} = strjoin(name, '');
%                 var_names_raw_data{end+1} = strjoin(new_column_x, '');
%                 var_names_raw_data{end+1} = strjoin(new_column_y, '');
%             end
%         end
%         i = i + 1;
%     end
% end
% 
% raw_table = readtable([paths.csvFolder filesep 'data_csv']);
% raw_table.Properties.VariableNames = var_names_raw_data; % Holds all raw data.
% smoothed_matrix = nan(height(raw_table), width(raw_table));
% excluded_matrix = nan(height(raw_table), width(raw_table)); % Will hold raw data without data points which fall in segement boundaries.
% smoothed_excluded_matrix = nan(height(raw_table), width(raw_table)); % Will hold data excluding segment boundaries fit with a smoothing spline.
% if isfile([paths.csvFolder filesep 'smoothed_table.csv'])
%     smoothed_table = readtable([paths.csvFolder filesep 'smoothed_table.csv']);
% end
% if isfile([paths.csvFolder filesep 'smoothed_excluded_table.csv'])
%     smoothed_excluded_table = readtable([paths.csvFolder filesep 'smoothed_excluded_table.csv']);
% end
disp('Done');

%% Align profiles
disp('------------Aligning Profiles by Amputation Plane------------')
f = figure;
g = figure;
for i = 1:width(analysis_mat)
    name = ['fish' num2str(analysis_mat(i).fish) '_ray' num2str(analysis_mat(i).ray) '_' num2str(analysis_mat(i).hpa) 'hpa'];
    disp(name);
    
    if analysis_mat(i).hpa == 0
        reference_amp_x = analysis_mat(i).amp_plane(1, 1);
    end
    
    amp_x = analysis_mat(i).amp_plane(1, 1);
    shift = amp_x - reference_amp_x;
    analysis_mat(i).x_pixels_shifted = analysis_mat(i).x_pixels + shift;
    analysis_mat(i).x_microns_shifted = analysis_mat(i).x_pixels_shifted * x_conversion;
    figure(f); plot(analysis_mat(i).x_microns_shifted, analysis_mat(i).raw_profile); hold on;
    title('aligned');
    figure(g); plot(analysis_mat(i).x_microns, analysis_mat(i).raw_profile); hold on;
    title('raw)');
end
    
%% Find and Exclude Segment Boundaries
disp('------------Identifying and Excluding Segment Boundaries------------');
j = 1;
for k = 1:width(analysis_mat)
    curr_im = analysis_mat(k);
    name = ['fish' num2str(curr_im.fish) '_ray' num2str(curr_im.ray) '_' num2str(curr_im.hpa) 'hpa'];
    disp(name);
    % if startsWith(name, 'fish1_ray6_3') || startsWith(name, 'fish3_ray6_1')
    %    continue;
    % else
    % x_name = strcat(name, '_x');
    % y_name = strcat(name, '_y');
    % raw_x = raw_table{:, x_name};
    % raw_y = raw_table{:, y_name};
    % h = height(raw_y(~isnan(raw_y))); % Height of the signal excluding all values which equal NaN.
    % raw_x = raw_x(1:h);
    % raw_y = raw_y(1:h);
    % smooth_raw_y = feval(fit_spline(raw_x, raw_y), raw_x);

    window = 100;
    if curr_im.hpa == 0 % Only finds boundaries on profiles from 0hpa and propagates to all the subsequent timepoints.
        boundaries = findSegmentBoundaries(curr_im.x_pixels, curr_im.raw_profile, window);
    end
    
    new_y = curr_im.raw_profile(1:boundaries(1));
    for i = 1:height(boundaries)/2
        x1 = boundaries(i*2, :);
        if i == height(boundaries)/2
            x2 = height(curr_im.x_pixels);
        else
            x2 = boundaries((i*2)+1, :);
        end
        new_y = [new_y; NaN(x1-height(new_y)-1, 1); curr_im.raw_profile(x1:x2)];
    end
    new_x = curr_im.x_pixels;
    smooth_y = feval(fit_spline(new_x, new_y(1:height(new_x))), new_x);

    curr_im.excluded_x_pixels = new_x;
    curr_im.excluded_x_microns = curr_im.excluded_x_pixels * x_conversion;
    curr_im.excluded_y = new_y;
    curr_im.smooth_excluded_y = smooth_y;
    % smoothed_matrix(1:height(curr_im.x_pixels), j*2-1) = curr_im.x_pixels;
    % smoothed_matrix(1:height(curr_im.raw_profile), j*2) = smooth_curr_im.raw_profile;
    % excluded_matrix(1:height(new_x), j*2-1) = new_x;
    % excluded_matrix(1:height(new_y), j*2) = new_y;
    % smoothed_excluded_matrix(1:height(new_x), j*2-1) = new_x;
    % smoothed_excluded_matrix(1:height(smooth_y), j*2) = smooth_y;

    f = figure;
    figure(f);
    plot(curr_im.x_pixels, curr_im.raw_profile);
    hold on
    plot(curr_im.excluded_x_pixels, curr_im.smooth_excluded_y);
    for bounds = boundaries(:, 1)%*x_conversion
        xline(bounds);
        hold on
    end
    hold off
    title(name, 'Interpreter', 'none');
    xlim([1, 1200]);
    ylim([1, 255]);
    % saveas(f, strjoin([paths.boundariesFolder filesep name '_boundaries.png'], ''));
    % close;

    j = j+1;
    % end
end

% smoothed_table = array2table(smoothed_matrix);
% excluded_table = array2table(excluded_matrix);
% smoothed_excluded_table = array2table(smoothed_excluded_matrix);
% smoothed_table.Properties.VariableNames = var_names_raw_data;
% excluded_table.Properties.VariableNames = var_names_raw_data;
% smoothed_excluded_table.Properties.VariableNames = var_names_raw_data;
% writetable(smoothed_table, [paths.csvFolder filesep 'smoothed_table.csv']);
% writetable(excluded_table, [paths.csvFolder filesep 'excluded_table.csv']);
% writetable(smoothed_excluded_table, [paths.csvFolder filesep 'smoothed_excluded_table.csv']);
disp('Done');

%% Plot Smoothened Ray Profiles
disp('------------Plotting Smoothened Profiles------------');
set(groot,'defaultLineLineWidth',2.0)
input_max = 800;
output_max = 255;
newcolors = brewermap(width(timepoints), '+Reds');

smoothed_table = readtable([paths.csvFolder filesep 'smoothed_table.csv']);
smoothed_excluded_table = readtable([paths.csvFolder filesep 'smoothed_excluded_table.csv']);

i = 1;
for fish = fish_nums
    for ray = ray_nums
        for i = 1:num_bifurcates(i)
            name = ['fish' fish '_ray' ray '_' num2str(i)];
            smoothed_savename = strjoin([paths.plotsFolder filesep 'fish' fish '_ray' ray '_' num2str(i) 'smoothed'], '');
            smoothed_video = VideoWriter(smoothed_savename,'MPEG-4'); %open video file
            smoothed_video.FrameRate = 1;
            open(smoothed_video)
            g = figure; % Figure corresponding to grouped plots
            figure(g);
            colororder(newcolors);
            
            smoothed_excluded_savename = strjoin([paths.plotsFolder filesep 'fish' fish '_ray' ray '_' num2str(i) 'smoothed_excluded'], '');
            smoothed_excluded_video = VideoWriter(smoothed_excluded_savename,'MPEG-4'); %open video file
            smoothed_excluded_video.FrameRate = 1;
            open(smoothed_excluded_video)
            f = figure; % Figure corresponding to grouped plots
            figure(f);
            colororder(newcolors);
            for time = timepoints
                x_name = strjoin([name '_' time 'hpa' '_x'], '');
                y_name = strjoin([name '_' time 'hpa' '_y'], '');
                smoothed_x = smoothed_table{:, x_name};
                smoothed_y = smoothed_table{:, y_name};                
                smoothed_excluded_x = smoothed_excluded_table{:, x_name};
                smoothed_excluded_y = smoothed_excluded_table{:, y_name};

                figure(g);
                plot(smoothed_x*x_conversion, smoothed_y);
                xlim([0 1000]);
                ylim([0 255]);
                title(strjoin(name, ''), 'Interpreter', 'none');
                xlabel('Distance from the amputation plane (microns)');
                ylabel('Pixel intensity (AU)');
                legend('0hpa', '4hpa', '8hpa', '12hpa', '16hpa', '20hpa', '24hpa'); %legend('0hpa', '12hpa', '24hpa', '48hpa', '60hpa', '72hpa'); %
                hold on;
                frame = getframe(g); %get frame
                writeVideo(smoothed_video, frame);
                
                figure(f);
                plot(smoothed_excluded_x*x_conversion, smoothed_excluded_y);
                xlim([0 1000]);
                ylim([0 255]);
                title(strjoin(name, ''), 'Interpreter', 'none');
                xlabel('Distance from the amputation plane (microns)');
                ylabel('Pixel intensity (AU)');
                legend('0hpa', '4hpa', '8hpa', '12hpa', '16hpa', '20hpa', '24hpa'); %legend('0hpa', '12hpa', '24hpa', '48hpa', '60hpa', '72hpa'); %
                hold on;
                frame = getframe(f); %get frame
                writeVideo(smoothed_excluded_video, frame);
            end
            close(g);
            close(smoothed_video);
            close(f);
            close(smoothed_excluded_video);
        end
        i = i + 1;
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