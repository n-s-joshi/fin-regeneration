%% preparation
%%%%
% this pipeline is compatible with MATLAB2021b
clear;clc;close all;
addpath(genpath('/Volumes/NSJ_Data_I/scripts/tgmm_alvin/segmentation_package/'))
addpath(genpath('./functions/'))


%% Initialize
dataFolder = 'F:\caudal_fin\04092024_osx-H2A-mCh_runx2a-GFP';
dataFolder = '/Volumes/NSJ_Data_I/caudal_fin/05_03_2025_Runx2-GFP';

dataFolder = [dataFolder filesep];
dataFolder_split = split(dataFolder,filesep);

% winFolder = 'C:\Users\Ashley\Desktop'; % directory to put data on workstation for TGMM run
winFolder = char(join(dataFolder_split(1:end-2),filesep));

TGMMtemplateFolder = 'F:\scripts\tgmm_alvin';
TGMMtemplateFolder = '/Volumes/NSJ_Data_I/scripts/tgmm_alvin';


TGMMtemplateFolder = [TGMMtemplateFolder filesep];

TGMMbuildFolder = 'F:\scripts\tgmm_alvin\TGMM_Supplementary_Software_1_0\build';
TGMMbuildFolder = '/Volumes/NSJ_Data_I/scripts/tgmm_alvin/TGMM_Supplementary_Software_1_0/build';

readTGMMFolder = ['./readTGMM_XMLoutput/'];

% bckgrd = 10;
% tau = 8;

ch_of_interest = [1]; % channels to quantify
ch_merge = [1,2]; % can have a merged channel if needed, if not leave as empty, the merged channel becomes channel 0. edited 070324
ch_ref = [1]; % reference channel for preprocessing
ch_nuc = [];
ch_KTR = []; % put [] if no KTR channel
ch_GEM = []; % put [] if no GEM channel
ch_BF = [2]; % put [] if no Bright Field channel
osteoblast = 0;
geometry = 1; % whether to calulate crossection area
rotate_YX = 1; % edited 092324 whether to rotate in YX plane
rot_or_not_post_tgmm = 0; % edited 092724 whether to rotate in YX plane after tgmm

% pluckTime={10 11 13 14 16 17; ...
%     datetime('2023-05-16 09:00:00') datetime('2023-05-16 09:00:00') datetime('2023-05-17 10:30:00')...
%     datetime('2023-05-17 10:30:00') datetime('2023-05-18 10:30:00') datetime('2023-05-18 10:30:00')}; %time of scale plucking
pluckTime= datetime('2024-11-10 11:00:00');

% ch_nuc = [2];
% ch_KTR = [1]; % put [] if no KTR channel
% ch_GEM = [4]; % put [] if no GEM channel
%% copy file (use this when there's two rays in one movie
%%
paths=[];
% paths.masterFolder='/Volumes/AshleyData2/28Mar22
% _H2A_ERK_Gem_fins_MEKtitrations/'; %folder where data is stored
paths.masterFolder=dataFolder; %folder where data is stored

paths.inFolder=  [paths.masterFolder 'raw/']; %input folder


fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*']);

% generate no_copy_list
noduplicate = ['4-2-72 4-6-72 2-2-108 4-2-108 4-6-108 2-2-168' ...
    ' 4-2-168 4-6-168 2-2-264 2-6-264 4-2-264 4-6-264 6-2-264'];

% noduplicate = strrep(noduplicate,'*','.'); % wild card for regexp

split_by_space = strsplit(noduplicate,' ');
no_copy_list = {};
for j = 1:numel(split_by_space)
    frt = strsplit(split_by_space{j},'-');
    no_copy_list = [no_copy_list,['fish' frt{1} '_ray' frt{2} '_' frt{3} 'hpa']];
end

for i=1:numel(st_dir)
    if contains(st_dir(i).name,'_ray2_')&&(~exist([paths.inFolder filesep strrep(st_dir(i).name,'ray2','ray3')],'dir'))
        disp(st_dir(i).name)      
        % no_copy_list = {'fish1_ray2_312hpa','fish2_ray2_312hpa','fish6_ray2_312hpa','fish1_ray2_384hpa',...
%             'fish2_ray2_384pa','fish2_ray2_552hpa'};
%             regexptranslate('wildcard','fish*_ray6_360hpa*')};
        if ~any(cell2mat(regexp(st_dir(i).name,no_copy_list)))
%         if ~ismember(st_dir(i).name,no_copy_list)
            copyfile([paths.inFolder filesep st_dir(i).name], ...
                [paths.inFolder filesep strrep(st_dir(i).name,'ray2','ray3')])
        end
    end
    if contains(st_dir(i).name,'_ray6_')&&(~exist([paths.inFolder filesep strrep(st_dir(i).name,'ray6','ray7')],'dir'))
        disp(st_dir(i).name)      
        % no_copy_list = {'fish3_ray6_312hpa'};
        if ~any(cell2mat(regexp(st_dir(i).name,no_copy_list)))
%         if ~ismember(st_dir(i).name,no_copy_list)
            copyfile([paths.inFolder filesep st_dir(i).name], ...
                [paths.inFolder filesep strrep(st_dir(i).name,'ray6','ray7')])
        end
    end
end

%% renaming file in each subfolders (use this when there's two rays in one movie
%
paths=[];
% paths.masterFolder='/Volumes/AshleyData2/28Mar22
% _H2A_ERK_Gem_fins_MEKtitrations/'; %folder where data is stored
paths.masterFolder=dataFolder; %folder where data is stored

paths.inFolder=  [paths.masterFolder 'raw/']; %input folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*']);

for i=1:numel(st_dir)
    if contains(st_dir(i).name,["_ray3_"])
        frt_dir = dir([paths.inFolder filesep st_dir(i).name filesep '**']); %fish ray time
        for j = 1:numel(frt_dir)
            if contains(frt_dir(j).name,["_ray2_"])
                [frt_dir(j).folder filesep frt_dir(j).name]
                movefile([frt_dir(j).folder filesep frt_dir(j).name], ...
                    [frt_dir(j).folder filesep strrep(frt_dir(j).name,'ray2','ray3')])
            end
        end
    end
    if contains(st_dir(i).name,["_ray7_"])
        frt_dir = dir([paths.inFolder filesep st_dir(i).name filesep '**']); %fish ray time
        for j = 1:numel(frt_dir)
            if contains(frt_dir(j).name,["_ray6_"])
                [frt_dir(j).folder filesep frt_dir(j).name]
                movefile([frt_dir(j).folder filesep frt_dir(j).name], ...
                    [frt_dir(j).folder filesep strrep(frt_dir(j).name,'ray6','ray7')])
            end
        end
    end
end

%% fix stitching define groups
to_fix = ['1-2-240 1-2-288 1-2-336 1-2-360 1-2-408 1-2-504 ' ...
    '1-3-240 1-3-288 1-3-336 1-3-360 1-3-408 1-3-504 1-6-192 1-6-240 1-6-288 1-6-312 1-6-480 ' ...
    '1-7-192 1-7-240 1-7-288 1-7-312 1-7-480 ' ...
    '2-6-240 2-6-480 2-7-216 2-7-408 ' ...
    '3-3-168 3-3-216 3-6-240 3-6-408 3-6-432 3-6-480 3-6-504 3-7-240 3-7-408 3-7-432 3-7-480 3-7-504 ' ...
    '4-6-384 4-7-408 5-2-312 5-3-312 5-2-432 5-3-432 5-6-240 5-6-288 5-6-360 5-7-240 5-7-288 5-7-360 ' ...
    '6-3-408 6-6-288 6-6-408 6-7-408'];
split_by_space = strsplit(to_fix,' ');
fix_list = {};
for j = 1:numel(split_by_space)
    frt = strsplit(split_by_space{j},'-');
    fix_list = [fix_list,['fish' frt{1} '_ray' frt{2} '_' frt{3} 'hpa']];
end
%% Stitching
% Stitch scales part
% The regenerating scale is usually image in 4-9 position that require
% need to be stitched in 3D. This section stitches the positions in a
% single stack.

close all;


paths=[];
% paths.masterFolder='/Volumes/AshleyData2/28Mar22
% _H2A_ERK_Gem_fins_MEKtitrations/'; %folder where data is stored
paths.masterFolder=dataFolder; %folder where data is stored

paths.inFolder=  [paths.masterFolder 'raw/']; %input folder 
paths.outFolder= [paths.masterFolder 'stitched/']; %output folder
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder

opts=[];
opts.stellaris = 1;
%opts.refCh=[ch_KTR ch_KTR ch_KTR]; %reference channel to use for stitching, can now be a list corresponding to each position, edited 091323
opts.refCh=[1]; %reference channel to use for stitching, can now be a list corresponding to each position, edited 091323
opts.targetVox=[0.6055,0.6055,0.6055]; %desidered voxel size (resizing may be applied)
opts.cropImage=0; %crop both stacks during stitching (0:no, 1:yes)
opts.verbose=0; %display images while stitching (0:no, 1:yes)
opts.inFolders=1; %files from the same experiment are contained in root folders with the same name (0:no, 1:yes)
opts.allSubs=0; % consider all position for stitching (0:no, 1:yes)
opts.timeStep = 1;
opts.resizeFactors=[0.2 0.2]; %image resizing factor used for stitching
opts.zBest=[]; % z to use for stitching ([position z]; [] for automatic choice)
opts.meta_pos=[1]; % whether or not to use metadata position to guide xcorr, can now be a list corresponding to each position, edited 020524
opts.laplacian_corr = 1; % use laplacian of the corr matrix for finding the stitching coordinates
opts.force_leftwards = [1]; % force stitching to go from right to left % edited 111723, can now be a list, edited 020914
opts.skip_save_error = 0; % whether to skip error and continue running and save errors in a log file

mkdir(paths.objFolder);
mkdir(paths.outFolder);

fish='';
scale=''; 
hpp='';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*_*' num2str(hpp) 'hpa']);

tic
for i=1:numel(st_dir)

   disp('------------');
   name=st_dir(i).name;    
    
   display([st_dir(i).name]);
          
            folder=[paths.inFolder st_dir(i).name '/'];
            basename=st_dir(i).name;
            % edited 120222, only run on data thats not been processed
            if exist([paths.objFolder basename '.mat'],'file')
                continue
            end
            % fish_to_skip = 'fish3_ray20_48hpa/';
            % if endsWith(folder, 'fish3_ray20_48hpa/')
            %     continue
            % end
            % if endsWith(folder, 'fish5_ray20_24hpa/')
            %     continue
            % end
            % if endsWith(folder, 'fish4_ray20_72hpa/')
            %     continue
            % end
            if(~exist(folder,'dir'))
                continue
            end

                if(exist([folder 'MetaData/' basename '_zBest.txt'],'file'))
                    opts.zBest=textread([folder 'MetaData/' basename '_zBest.txt']);
                elseif(exist([folder 'MetaData/zBest.txt'],'file'))
                    opts.zBest=textread([folder 'MetaData/zBest.txt']);
                end
                paths.rayFolder = [paths.outFolder filesep 'show_each_step' filesep basename];
                mkdir(paths.rayFolder);
                
                % set pluckTime for each ray, edited 052924
                if numel(pluckTime)>1
                    pluckTime_here = pluckTime{2,[pluckTime{1,:}]==sscanf(basename,'fish%d_ray*')};
                else
                    pluckTime_here = pluckTime;
                end

                if isfield(opts,'skip_save_error') && opts.skip_save_error
                    try % edited 120222
                        [myScale] = stitchScale(paths,basename,opts);
                        % write correct hpp
                        myScale.pluckTime=pluckTime_here;
                        myScale.hppTrue=myScale.metadata{1}.hTimeStamp-datenum(pluckTime_here)*24;
                        wpathmat=[paths.objFolder basename '.mat'];  
                        save(wpathmat,'myScale'); 
                    catch err
                        %open file
                        fid = fopen([paths.outFolder 'errorlog.txt'],'a+');
                        % write the error to file
                        % 0 line: fish ray hpa
                        fprintf(fid, [st_dir(i).name '\n']);
                        % first line: message
                        fprintf(fid,'%s\n',err.message);
                        % following lines: stack
                        for e=1:length(err.stack)
                        fprintf(fid,'%s in %s at line %i \n',err.stack(e).name,err.stack(e).file,err.stack(e).line);
                        end
                        fprintf(fid,'\n');
    
                        % close file
                        fclose(fid);
                    end
                else
                    [myScale] = stitchScale(paths,basename,opts);
                    % write correct hpp
                    myScale.pluckTime=pluckTime_here;
                    % myScale.hppTrue=myScale.metadata{1}.hTimeStamp-datenum(pluckTime_here)*24;
                    wpathmat=[paths.objFolder basename '.mat'];  
                    save(wpathmat,'myScale');
                end

%                 close all


end
toc

%% duplicate files
%%
paths=[];
% paths.masterFolder='/Volumes/AshleyData2/28Mar22
% _H2A_ERK_Gem_fins_MEKtitrations/'; %folder where data is stored
paths.masterFolder=dataFolder; %folder where data is stored

paths.inFolder=  [paths.masterFolder 'stitched/']; %input folder
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*']);

% generate no_copy_list
toduplicate = ['1-2-144 1-2-156 1-2-168'];
changerayto = '23';
% noduplicate = strrep(noduplicate,'*','.'); % wild card for regexp

split_by_space = strsplit(toduplicate,' ');
to_copy_list = {};
for j = 1:numel(split_by_space)
    frt = strsplit(split_by_space{j},'-');
    to_copy_list = [to_copy_list,['fish' frt{1} '_ray' frt{2} '_' frt{3} 'hpa']];
end

for i=1:numel(st_dir)
    if contains(st_dir(i).name,'_ray2_')&&(~exist([paths.inFolder filesep strrep(st_dir(i).name,'ray2',['ray' changerayto])],'file'))
        disp([st_dir(i).name '_stitched'])      
        % no_copy_list = {'fish1_ray2_312hpa','fish2_ray2_312hpa','fish6_ray2_312hpa','fish1_ray2_384hpa',...
%             'fish2_ray2_384pa','fish2_ray2_552hpa'};
%             regexptranslate('wildcard','fish*_ray6_360hpa*')};
        if any(cell2mat(regexp(st_dir(i).name,to_copy_list)))
            disp(['----copy----'])   
%         if ~ismember(st_dir(i).name,no_copy_list)
            copyfile([paths.inFolder filesep st_dir(i).name], ...
                [paths.inFolder filesep strrep(st_dir(i).name,'ray2',['ray' changerayto])])
        end
    end
end


% duplicate objects
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
for i=1:numel(st_dir)

    disp([st_dir(i).name '_objects'])  
    basename=st_dir(i).name;    
    display(basename);
    if any(cell2mat(regexp(basename,to_copy_list)))
        disp(['----copy----'])   
        load([st_dir(i).folder filesep basename]);
        myScale.basename=strrep(myScale.basename,'ray2',['ray' changerayto]);
        basename_out = strrep(basename,'ray2',['ray' changerayto]);
        wpathmat=[paths.objFolder basename_out];  
        save(wpathmat,'myScale');
    end
end



%% Change plucktime - only use this if you forgot to change plucktime in the intialization step and want to correct them
paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
 
for i=1:numel(st_dir)

    disp('------------');
    basename=st_dir(i).name;    
    display(basename);
    load([st_dir(i).folder filesep basename]);

    myScale.pluckTime=pluckTime;
    myScale.hppTrue=myScale.metadata{1}.hTimeStamp-datenum(pluckTime)*24;
    wpathmat=[paths.objFolder basename];  
    save(wpathmat,'myScale');
end


%% Add time, scale and fish
% This section reads the time-point filename and stores information on fish
% number, scale number and time labels.

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
 
for i=1:numel(st_dir)

   disp('------------');
   name=st_dir(i).name;    
    
   display([st_dir(i).name]);
   load([st_dir(i).folder filesep st_dir(i).name]);
    
   name=myScale.basename;
   
   delimiters=[regexp(name,'_') numel(name)+1];
   % time - hpp
   poshpp=regexp(name,'hpa');
   if(~isempty(poshpp))
        postime=delimiters(find(delimiters<poshpp,1,'last'))+1;
        myScale.hpp =str2num(name(postime:poshpp-1));
   else
        myScale.hpp =[]; 
   end
   
   % time - dpp
   posdpp=regexp(name,'hpa');
   if(~isempty(posdpp))
        postime=delimiters(find(delimiters<posdpp,1,'last'))+1;
        myScale.dpp =str2num(name(postime:posdpp-1));
        myScale.hpp = myScale.dpp*24;
   else
        myScale.dpp =[]; 
   end

   % time - t
   post=regexp(name,'_t');
   if(~isempty(post))
        posEndtime=delimiters(find(delimiters>post,1,'first'))-1;
        myScale.t =str2num(name(post+2:posEndtime));
   else
        myScale.t=[]; 
   end
   
   % Mark and Find
   posMF=regexp(name,'Mark_and_Find_');
   if(~isempty(posMF))
       posEndMF=delimiters(find(delimiters>posMF,4,'first'))-1;
       posEndMF=posEndMF(end);
       myScale.MF=str2num(name(posMF+14:posEndMF));
   else
       myScale.MF=[];
   end

   % fish
   posfish=regexp(name,'fish');
   posnumber=delimiters(find(delimiters>posfish,1,'first'))-1;
   fish =name((posfish+4):posnumber);

   % scale
   posray=regexp(name,'ray');
   posnumber=delimiters(find(delimiters>posray,1,'first'))-1;
   ray =name((posray+3):posnumber);

   myScale.fish = str2num(fish);
   myScale.scale= str2num(scale);
   myScale.ray= str2num(ray);


   save([st_dir(i).folder filesep st_dir(i).name],'myScale');
   
end
%% merge channels if needed
if numel(ch_merge)>1
disp(["---------- Merge reference channels ----------"]);

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=     [paths.masterFolder 'stitched/']; % input folder stack to equalize
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder




fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);


for i=1:numel(st_dir)
    display([st_dir(i).name]);
    load([st_dir(i).folder filesep st_dir(i).name]);
    
    if myScale.hpp==672 % out of ram
        % continue
    end


    basename=myScale.basename;
    if exist([paths.inFolder basename '_ch' num2str(0) 'maxproj.tif'],'file') % skip exist
        % continue
    end

    
    stack_4d = [];
    for ch=ch_merge
      
        spath=[paths.inFolder basename '_ch' num2str(ch) '.tif'];
        
        %using loadtiff
        stackHere=loadtiff(spath);
        
        stack_4d = cat(4,stack_4d,stackHere);
    end
    
    stack_max4 = max(stack_4d,[],4);
    
    %write channel
    wpath=[paths.inFolder basename '_ch' num2str(0) '.tif']; 
    
    options=[];
    options.overwrite='true';
    options.compress='lzw';
    saveastiff(stack_max4, wpath, options);

    % write max proj
    stack_maxproj = max(stack_max4,[],3);
    wpath_proj=[paths.inFolder basename '_ch' num2str(0) 'maxproj.tif']; 
    
    options=[];
    options.overwrite='true';
    options.compress='lzw';
    saveastiff(stack_maxproj, wpath_proj, options);

end


end


%% Cleans up  scales - Create ROI file - USER INPUT REQUIRED (SEE BELOW)
% Scale stack include neighboring scales portions. Manual data curation is
% required to eliminate neighbouring scales. In this step, the user selects the
% region where the  selected scale is located in three different views. 
%
% This section requires user input. When data sample is tested, this
% section is skipped as we already provide parameters.
% If users want to chose parameters, please set useSavedParams= false

useSavedParams = false;

refCh=[1]; % reference channels to use for ROI selections, can put 2 channels and will show overlap with intensity adapted for each channel separately, -1 for merge all channels
chToClean=[];  % reference channels to use for manual cleaning step. Need to be empty here
roi=[]; % roi to use (stack; [] to ask user input)

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=    [paths.masterFolder 'stitched/']; % input folder for images to cure
paths.outFolder=   [paths.masterFolder 'cleaned/']; % output folder for ROIs
%paths.outFolderSaved=   [paths.masterolder 'cleanedSaved/']; % output folder for ROIs
paths.fastFolder = [paths.outFolder filesep 'ROIproj/'];
paths.objFolder=   [paths.masterFolder 'objects/']; %objects folder

if(~useSavedParams)
    
    fish='3'; % put the number of the fish/scale/hpp to process or '' to choose them all
    scale='20'; 
    hpp='120';
    st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' 'ch' num2str(refCh(1)) '.tif']);

    % sort st_dir; edited 121222
    [~,sort_idx] = sort_nat({st_dir.name});
    st_dir = st_dir(sort_idx);
    opts=[];
    opts.refCh = refCh;
    opts.type='manual'; 
    opts.cleanupPlanes = 0;
    opts.adj = [0 1]; % the user can change image display adjustment ( 0 < [min max] < 1)
    opts.BF = ch_BF; % whether to use bright field for cropping
    opts.fast = 1; % whether to use the faster code (selectROIstack_fast.m)
    opts.roiTight = []; % whether to define additional roi (usuallyYX) using other channels, put empty if none
    opts.auto_crop_z = 1; % whether to only draw in xy and automatically crop in z based on a threshold that you can tweak
    
    for i=1:numel(st_dir)
            % edited 121222, only run on data thats not been processed
            if exist([paths.fastFolder st_dir(i).name(1:end-8) '_ROIprojYX.tif'],'file')
                % continue
            end
        disp([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
        load([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
        [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts);
        close all

        % save myScale
%         myScale.pluckTime = pluckTime;
%         wpathmat=[paths.objFolder myScale.basename '.mat'];  
%         save(wpathmat,'myScale');

    end
else
    
%    mkdir(paths.outFolder);
%    copyfile([paths.outFolderSaved '*ROI.tif'],[paths.outFolder])
%    copyfile([paths.outFolderSaved '*ch' num2str(refCh) '*.tif'],[paths.outFolder])
%    disp('Saved parameters are used, please continue to next section'); 

end

%% start timer

tStart = tic;
%% Convert fast ROIproj into roi in 3D
disp('---------Starting converting rois----------');

refCh=ch_ref; % reference channel (not used since you are providing ROIs)
chToClean=ch_BF; % channels to apply cleaning to
roi=[];

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder  =  [paths.masterFolder 'stitched/']; % input folder stacks to rotate
paths.roiFolder  =  [paths.masterFolder 'cleaned/']; % folder provided ROIs
paths.outFolder =  [paths.masterFolder 'cleaned/']; % output folder
paths.objFolder =  [paths.masterFolder 'objects/']; % objects folder 
paths.fastFolder = [paths.roiFolder filesep 'ROIproj/'];

fish='';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

opts=[];
opts.type='manual';
% opts.auto_crop_z = 1; % whether to only draw in xy and automatically crop in z based on refCh
opts.refCh = refCh;

for i=1:numel(st_dir)
    if exist([paths.roiFolder st_dir(i).name(1:end-4) '_ROI.tif'],'file')
         % continue
    end
    % edited 120922 process fast projROI into 3d roi
    % edited 052523 to enable multiple roiTight
    roiYX_dir = dir([paths.fastFolder filesep st_dir(i).name(1:end-4) '_ROIprojYX*.tif']);
    for j = 1:numel(roiYX_dir)
        name_split_cell = strsplit(roiYX_dir(j).name,'_');
        suffix = strjoin(name_split_cell(4:end),'_');
        suffix_save = strjoin(name_split_cell(5:end),'_');
%     if exist([paths.fastFolder filesep st_dir(i).name(1:end-4) '_ROIprojYX.tif'],"file")
        options=[];
        options.overwrite='true';
        options.compress='lzw';
        disp([st_dir(i).name]);
        load([paths.objFolder filesep st_dir(i).name]);
        % run
        roi = convertROIprojto3d(myScale,paths,suffix,opts);
        % save roi image
        if isempty(suffix_save)
            wpath=[paths.roiFolder myScale.basename '_ROI.tif']; 
        else
            wpath=[paths.roiFolder myScale.basename '_ROI_' suffix_save]; 
        end
        myScale.cleanedScale.roipath=wpath;
        saveastiff(im2uint8(roi), wpath, options);
        % save myScale
        wpathmat=[paths.objFolder myScale.basename '.mat'];  
        save(wpathmat,'myScale');

    end
end
%% Apply ROI to clean (v.) images 
% Apply ROI3 to rotated images to get cleaned images

disp('---------Starting cleaning images----------');

refCh=ch_ref; % reference channel (not used since you are providing ROIs)
chToClean=ch_BF; % channels to apply cleaning to
roi=[];

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder  =  [paths.masterFolder 'stitched/']; % input folder stacks to rotate
paths.roiFolder  =  [paths.masterFolder 'cleaned/']; % folder provided ROIs
paths.outFolder =  [paths.masterFolder 'cleaned/']; % output folder
paths.objFolder =  [paths.masterFolder 'objects/']; % objects folder 
paths.fastFolder = [paths.roiFolder filesep 'ROIproj/'];

fish='';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir(strcat(paths.roiFolder,['fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' ],'ROI.tif'));

opts=[];
opts.type='manual';
 
for i=1:numel(st_dir)
    if exist([paths.outFolder st_dir(i).name(1:end-8) '_ch' num2str(chToClean(1)) '.tif'],'file')
        % continue
    end
    disp([st_dir(i).name(1:end-8) '.mat']);
    load([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
    
    % edited 120922 process fast projROI into 3d roi
%     if exist([paths.fastFolder filesep myScale.basename '_ROIprojYX.tif'],"file")
%         options=[];
%         options.overwrite='true';
%         options.compress='lzw';
%         roi = convertROIprojto3d(myScale,paths);
%         wpath=[paths.roiFolder myScale.basename '_ROI'  '.tif']; 
%         myScale.cleanedScale.roipath=wpath;
%         saveastiff(im2uint8(roi), wpath, options);
%     end

    roi=loadtiff([paths.roiFolder filesep st_dir(i).name]);
    [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts);

    % record cropping rectangle (3d)
    roi_ind = find(logical(roi));
    [roiIdx1,roiIdx2,roiIdx3] = ind2sub(size(roi),roi_ind);
    crop_rect = [max(min(roiIdx1)-10,1),min(max(roiIdx1)+10,size(roi,1));...
        max(min(roiIdx2)-10,1),min(max(roiIdx2)+10,size(roi,2));...
        max(min(roiIdx3)-2,1),min(max(roiIdx3)+2,size(roi,3))];
    myScale.cleanedScale.crop_rect = crop_rect;
   %save mat file
   wpathmat=[paths.objFolder myScale.basename '.mat'];  
   save(wpathmat,'myScale');
end

%% one-step crop
disp('----------Starting cropping image')

chToCrop=ch_BF; % channels to apply cleaning to

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder  =  [paths.masterFolder 'cleaned/']; % input folder stacks to rotate
paths.outFolder =  [paths.masterFolder 'cleaned/cropped/']; % output folder
paths.objFolder =  [paths.masterFolder 'objects/']; % objects folder 
paths.fieldName=  'cleanedScale'; %struct where re-rotation parameters are read

fish='';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir(strcat(paths.inFolder,['fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' ],'ROI.tif'));

opts=[];
opts.type='manual';
 
for i=1:numel(st_dir)
    if exist([paths.outFolder st_dir(i).name(1:end-8) '_ch' num2str(chToCrop(1)) '.tif'],'file')
%         continue
    end
    disp([st_dir(i).name(1:end-8) '.mat']);
    load([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
    
    [myScale] = crop_onestep(myScale,chToCrop,paths);
end

%% Resize cleaned and cropped scales 
% This section resizes cleaned images to fasten calculation of optimal
% re-rotation angles.
disp('----------Resize images for rotation---------')

chToResize = ch_BF; %reference channel for angles calculation

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder= [paths.masterFolder 'cleaned/cropped/']; % input folder cleaned images
paths.outFolder= [paths.masterFolder 'cleaned/resized4x/']; % output folder resized images
paths.objFolder= [paths.masterFolder 'objects/']; % objects folder
paths.fieldName=  'resizeCropped'; % struct field resizing parameters

fish='';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' 'ch' num2str(chToResize) '.tif']);

opts=[];
opts.resizeFactor = 0.25; % resizing factor 

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

mkdir(paths.outFolder);

for i=1:numel(st_dir)
    if exist([paths.outFolder st_dir(i).name(1:end-8) '_ch' num2str(chToResize(1)) '.tif'],'file')
%         continue
    end
    display([st_dir(i).name]);
    load([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
 
    for ch=1:numel(chToResize)
        s=loadtiff([paths.inFolder filesep st_dir(i).name(1:end)]);
        startSize=size(s);
        s1 = imresize(s,opts.resizeFactor);
       
        xSize = size(s1,1);
        s1=permute(s1,[1 3 2]);
        s2 = imresize(s1,[xSize round(opts.resizeFactor.*startSize(3))]);
        s2=permute(s2,[1 3 2]);
        saveastiff(s2,[paths.outFolder filesep st_dir(i).name(1:end)],optsSave);
    end
    
    myScale.(paths.fieldName).startSize=startSize;
    myScale.(paths.fieldName).resizeFactor=opts.resizeFactor;
    save([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat'],'myScale');

end
disp('Done');
%% Re-rotate scales - Calculate rotation parameters - (actually theres only one rotation in this pipeline)
% The stack is re-rotated so that the cleaned scale is as parallel as possibile
% to the xy plane. Angles are automatically calculated and the user chooses
% cropping rectangles to contain only the rotated scales. 
%
% This section requires user input. When data sample is tested, this
% section is skipped as we already provide parameters.
% If users want to chose parameters, please set useSavedParams= false
disp(["----------Starting calculate rotation parameters----------"]);

useSavedParams= false;

if(~useSavedParams)
    refCh=ch_ref;
    chToRotate = ch_of_interest; %channels to rotate at this step ([] for none)
    if rotate_YX % edited 092324 whether to rotate in YX plane
        ths=[NaN NaN NaN];
    else
        ths=[]; %rotation angles ([th1 th2 th3]; [] to calculate angles automatically)
    end
    rects=zeros(4,4);
    % cropping rectangles in different projections ([positionX positionY width height]; [NaN NaN NaN
    % NaN] to ask user input; [0 0 0 0] to not crop)
    zPos=0; % planes to use ([zMin zMax]; 0 to keep them all)

    paths=[];
    paths.masterFolder=dataFolder; %folder where data is stored
    % paths.inFolder=  [paths.masterFolder 'cleaned/resized4x/']; %images
    % to use to calculate rotation % changed to line 812 on 4/8/25
    paths.inFolder=  [paths.masterFolder 'cleaned/'];
    paths.outFolder= [paths.masterFolder 'rotated/']; % output folder (not used here)
    paths.objFolder= [paths.masterFolder 'objects/']; % objects folder

    paths.fieldName=  'resizeCropped'; %struct field where parameters are read
    paths.rotFieldName = 'rotScale';    %struct field where rotation parameters are stored

    fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
    scale=''; 
    hpp='';

    st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
    % sort st_dir; edited 062123
    [~,sort_idx] = sort_nat({st_dir.name});
    st_dir = st_dir(sort_idx);

    resizeFactor = 0.25; % resizing factor

    % calculate rotations
    for i=1:numel(st_dir)
        if exist([paths.outFolder st_dir(i).name(1:end-4) '_ch' num2str(ch_of_interest(1)) '.tif'],'file')
%             continue
        end
        if ~exist([paths.inFolder st_dir(i).name(1:end-4) '_ch' num2str(ch_ref) '.tif'],'file')
            % continue
        end
        display([st_dir(i).name]);
        load([paths.objFolder filesep st_dir(i).name]);

        myScale.(paths.fieldName).resizeFactor = resizeFactor;
        [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths);
        myScale.(paths.rotFieldName).rotAngles = myScale.(paths.fieldName).rotAngles;
        myScale.(paths.rotFieldName).zPos = round(myScale.(paths.fieldName).zPos./(myScale.(paths.fieldName).resizeFactor));

        myScale.(paths.rotFieldName).cropRects = rects;
        save([paths.objFolder filesep st_dir(i).name],'myScale');
    end
else
    disp('Saved parameters are used, please continue to next section'); 
end

%% Re-rotate scales - (actually theres only one rotation step (z) in this pipeline)
disp(["----------Starting rotation and crop to remove excessive boundaries----------"]);
refCh=ch_ref;
chToRotate = ch_of_interest; %channel to be processed

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=  [paths.masterFolder 'cleaned/cropped/']; % input folder for images to re-rotate
paths.outFolder= [paths.masterFolder 'rotated/'];  % output folder for re-rotated images
paths.objFolder= [paths.masterFolder 'objects/'];  % objects folder
%paths.objFolderSaved= [paths.masterFolder 'objectsSaved/']; % provided objects folder
paths.fieldName=  'rotScale'; %struct where re-rotation parameters are read

useSavedParams = false; % use provided parameters (true: use provided parameters; false: use user parameters)

fish='3'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='20'; 
hpp='120';

st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
% sort st_dir; edited 062123
[~,sort_idx] = sort_nat({st_dir.name});
st_dir = st_dir(sort_idx);

% apply re-rotations
for i=1:numel(st_dir)
    if exist([paths.outFolder st_dir(i).name(1:end-4) '_ch' num2str(chToRotate(1)) '.tif'],'file')
            % continue
    end
    if ~exist([paths.inFolder st_dir(i).name(1:end-4) '_ch' num2str(chToRotate(1)) '.tif'],'file')
            % continue
    end
    disp([st_dir(i).name]);
    
    %load([paths.objFolderSaved filesep st_dir(i).name]);
    %myScaleSaved = myScale;
    load([paths.objFolder filesep st_dir(i).name]);
   
    if(useSavedParams)
        ths   = myScaleSaved.(paths.fieldName).rotAngles; 
        rects = myScaleSaved.(paths.fieldName).cropRects;
        zPos  = myScaleSaved.(paths.fieldName).zPos;
    else
        ths   = myScale.(paths.fieldName).rotAngles; 
        rects = myScale.(paths.fieldName).cropRects;
        zPos  = myScale.(paths.fieldName).zPos;
    end
    
    [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths);
    save([paths.objFolder filesep st_dir(i).name],'myScale');
end


%% Equalization nuclear signal for segmentation and tracking 
disp(["----------Start equalization----------"]);


chToEq=ch_of_interest; % channel to equalize

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=     [paths.masterFolder 'rotated/']; % input folder stack to equalize
paths.outFolder=    [paths.masterFolder 'equalized/']; % output folder equalized stacks
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

 % equalization
for i=1:numel(st_dir)
    display([st_dir(i).name]);
    if exist([paths.outFolder st_dir(i).name(1:end-4) '_ch' num2str(chToRotate(1)) '.tif'],'file')
          continue
    end
    load([st_dir(i).folder filesep st_dir(i).name]);
    myScale = equalizeScale(myScale,chToEq,paths);
end

%% These blocks only required for osteoblasts
if osteoblast

%%  Entire scale segmentation in 3D (not needed for fibroblasts)
disp(["----------Start entire tissue segmentation----------"]);

refCh=ch_ref; % reference channels for segmentation
paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=     [paths.masterFolder 'rotated/']; % input folder stack to segment %cleanedRerotFlipped
paths.outFolder= [paths.masterFolder 'segmented/']; %segmented 3D ROIs
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder

%As of 1June21: Use cleanedRerot Ch3 with threshold of 2.5

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
verbose =0; % display segmentation (0: no; 1: yes)

segOpts.segmThresh = 1; %Use 1.0 on 22Dec21

for i=1:numel(st_dir)
    disp('------------');
    display([st_dir(i).name(1:end-4) '.mat']);
    if exist([paths.outFolder st_dir(i).name(1:end-4) '_ch' num2str(ch_of_interest(1)) '_seg.tif'],'file')
          continue
    end
    load([paths.objFolder filesep st_dir(i).name(1:end-4) '.mat']);
    [myScale] = segmentMyScale(myScale,refCh,paths,segOpts,verbose);
end

%% Scale flattening - Equalized stacks (not needed for fibroblasts)
% The scale is flattened by bringing its segmented hyposquamal border at the same
% z-position. 
disp(["----------Start flattening equalized stacks----------"]);

toflattenCh=ch_of_interest; % channel to flatten

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder= [paths.masterFolder 'equalized/']; % input folder equalized stacks
paths.roiFolder= [paths.masterFolder 'segmented/']; % segmented ROI folder
paths.outFolder= [paths.masterFolder 'flatten_eq/']; % output folder flattened stacks
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder
paths.fitFieldName='findHypo'; % struct field flattening parameters

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

opts=[];
opts.method='fit'; % with 'fit' ROI is not flattened - launch giving toFlattenCh='ROI' to flatten ROI
 
for i=1:numel(st_dir)
    
    display([st_dir(i).name]);
    load([st_dir(i).folder filesep st_dir(i).name]);  
    [myScale] = flattenLayersScaleSegm(myScale,toflattenCh,paths,opts);
    
end
%% Scale flattening - Not equalized stacks (not needed for fibroblasts)
% The scale is flattened by bringing its segmented hyposquamal border at the same
% z-position. 
disp(["----------Start flattening un-equalized stacks----------"]);

toflattenCh=ch_of_interest; % channel to flatten

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=     [paths.masterFolder 'rotated/']; % input folder stacks
paths.roiFolder= [paths.masterFolder 'segmented/']; % segmented ROI folder
paths.outFolder= [paths.masterFolder 'flatten/']; % output folder flattened stacks
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder
paths.fitFieldName='findHypo';

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

opts=[];
opts.method='fit'; % with 'fit' ROI is not flattened - launch giving toFlattenCh='ROI' to flatten ROI
 
for i=1:numel(st_dir)
    
    display([st_dir(i).name]);
    load([st_dir(i).folder filesep st_dir(i).name]);  
    [myScale] = flattenLayersScaleSegm(myScale,toflattenCh,paths,opts);
    
end
%% Isolation hyposquamal layer - Equalized stacks (not needed for fibroblasts)
% The hyposquamal layer is computationally dissected . 
% To this end, the peak of the total plane intensity
% is calculated. The hyposquamal layer is taken from the start of the stack 
% to 15 planes past the z-peak of total plane intensity.
disp(["----------Start isolating hypo-layer for equalized stack----------"]);

refCh=ch_ref; % reference channel for calculation z-peak
todivideCh=ch_of_interest; % channels to divide

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=     [paths.masterFolder 'flatten_eq/']; % input folder stacks to divide
paths.refFolder=    [paths.masterFolder 'flatten_eq/']; % reference folder stack to calculate z-peak
paths.outFolder=    [paths.masterFolder 'divided_eq/']; % output folder for isolated hyposquamal layers
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folders
paths.fitFieldName='findHypo'; % struct field name for hyposquaml layer isolation parameters

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

divOpts=[];
divOpts.verbose=0;    % display z-intensities
divOpts.methodEpi=15; % z-limit hyposquamal layer - planes past the peak in total z-intensity

for i=1:numel(st_dir)
    
    display([st_dir(i).name]);
    if exist([paths.outFolder st_dir(i).name(1:end-4) '_ch' num2str(ch_of_interest(1)) '_sharedmaxproj.tif'],'file')
        continue
    end

    load([paths.objFolder st_dir(i).name]);

    % if(myScale.hpp<48)
    %     divOpts.methodHypo='peak';
    %     [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    % elseif(myScale.hpp<85)
    %     divOpts.methodHypo='peak';
    %     [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    % else
    %     divOpts.methodHypo='peak';
    %     [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    % end
    
    divOpts.methodHypo='peak';
    [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
end

%% Isolation hyposquamal layer - Non equalized stacks (not needed for fibroblasts)
% The hyposquamal layer is computationally dissected . 
% To this end, the peak of the total plane intensity
% is calculated. The hyposquamal layer is taken from the start of the stack 
% to 15 planes past the z-peak of total plane intensity.
disp(["----------Start isolating hypo-layer for un-equalized stack----------"]);

refCh=ch_ref;
todivideCh=ch_of_interest;

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.inFolder=     [paths.masterFolder 'flatten/'];  % input folder stacks to divide
paths.refFolder=    [paths.masterFolder 'flatten_eq/']; % reference folder stack to calculate z-peak
paths.outFolder=    [paths.masterFolder 'divided_refEq/']; % output folder for isolated hyposquamal layers
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folders

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

divOpts=[];
divOpts.verbose=0;
divOpts.methodEpi=15; % z-limit hyposquamal layer - planes past the peak in total z-intensity
  
for i=1:numel(st_dir)
    
    display([st_dir(i).name]);
    load([paths.objFolder st_dir(i).name]);

    % if(myScale.hpp<48)
    %     divOpts.methodHypo='peak';
    %     [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    % elseif(myScale.hpp<80)
    %     divOpts.methodHypo='peak';
    %     [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    % else
    %     divOpts.methodHypo='peak';
    %     [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    % end

    divOpts.methodHypo='peak';
    [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
   
end
%%
end

%% Get intensity profiles
% Written by NSJ on 03/11/2025
% Run to manually straighten and get the intensity profiles of images
disp('----------Begin Get Intensity Profiles----------');
paths = [];
paths.masterFolder = dataFolder;
paths.inFolder = [paths.masterFolder filesep 'rotated/'];
paths.outFolder = [paths.masterFolder filesep 'straightened/'];
mkdir(paths.outFolder);
paths.objFolder = [paths.masterFolder filesep 'objects/'];
chToUse = ch_of_interest;

fish = '3'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale = '2';
hpp = '60';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '_*' num2str(hpp) 'hpa*' '.mat']);
[~,sort_idx] = sort_nat({st_dir.name});
st_dir = st_dir(sort_idx);

for i=1:numel(st_dir)
    disp([st_dir(i).name]);
    wpathmat = [paths.objFolder filesep st_dir(i).name(1:end-4) '.mat'];
    load(wpathmat);
    if isfield(myScale, 'straighten_points')
        % continue
    end
    fluor_image_path = [paths.inFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(chToUse) 'maxproj.tif'];
    fluor_IM = imread(fluor_image_path);
    fluor_IM = widenImage(fluor_IM);
    points = getStraightenPoints(fluor_IM);
    % bf_image_path = [paths.inFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(ch_BF) 'maxproj.tif'];
    % bf_IM = imread(bf_image_path);
    % bf_IM = widenImage(bf_IM);
    myScale.straighten_points = points;
    w = 200;
    [fluor_IM2, w] = optimizeStraighten(fluor_IM, points', w);
    % bf_IM2 = straighten(bf_IM, points', w);
    close;
    imshow(permute(fluor_IM2,[2,1,3])./255); axis image off;
    save_path = [paths.outFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(chToUse) 'maxproj_straightened.tif'];
    saveas(gcf, save_path);
    close;
    % imshow(permute(bf_IM2,[2,1,3])./255); axis image off;
    % save_path = [paths.outFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(ch_BF) 'maxproj_straightened.tif'];
    % saveas(gcf, save_path);
    % close;
    profile = collapseMatrix(fluor_IM2);
    myScale.linescan = flip(profile);
    save(wpathmat, 'myScale');
end
disp('done');

%% Manually determine amputation plane
% Written by NSJ on 2/10/25
% Run to manually determine the location of the amputation plane in each
% ray
disp('----------Begin Manually Determine Amputation Plane----------');
paths = [];
paths.masterFolder = dataFolder;
paths.inFolder = [paths.masterFolder filesep 'straightened'];
paths.inFolder = [paths.masterFolder filesep 'rotated'];
paths.objFolder = [paths.masterFolder filesep 'objects'];
chToUse = ch_GEM;

fish = '9'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale = '6';
hpp = '';
st_dir=dir([paths.objFolder filesep 'fish'  num2str(fish) '_ray' num2str(scale) '_*' num2str(hpp) 'hpa*' '.mat']);
[~,sort_idx] = sort_nat({st_dir.name});
st_dir = st_dir(sort_idx);

for i=1:numel(st_dir)
    disp([st_dir(i).name]);
    wpathmat = [paths.objFolder filesep st_dir(i).name(1:end-4) '.mat'];
    load(wpathmat);
    if isfield(myScale, 'manualAmpPlane')
        % continue
    end
    file_end_name = 'maxproj_straightened.tif';
    file_end_name = 'maxproj.tif';
    proj = imread([paths.inFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(chToUse) file_end_name]);
    h = figure;
    h.Name = [st_dir(i).name(1:end-4) ' BF'];
    imagesc(proj); axis image off;
    h.WindowState = 'maximized';
    [x, y] = ginput(1);
    myScale.manualAmpPlane = [x, y];
    % myScale.straightenedXLength = width(proj); 
    save(wpathmat, 'myScale');
    close;
end

%% Find blastemal signal area
% Written by NSJ on 5/23/25
% Run to determine the area of a signal of interest in the blastema 
disp('----------Begin Find Blastemal Signal Area----------');
paths = [];
paths.masterFolder = dataFolder;
paths.inFolder = [paths.masterFolder filesep 'rotated'];
paths.objFolder = [paths.masterFolder filesep 'objects'];
chToUse = ch_ref;

fish = ''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale = '';
hpp = '';
st_dir=dir([paths.objFolder filesep 'fish'  num2str(fish) '_ray' num2str(scale) '_*' num2str(hpp) 'hpa*' '.mat']);
[~,sort_idx] = sort_nat({st_dir.name});
st_dir = st_dir(sort_idx);

for i=1:numel(st_dir)
    disp([st_dir(i).name]);
    wpathmat = [paths.objFolder filesep st_dir(i).name(1:end-4) '.mat'];
    load(wpathmat);
    if isfield(myScale, 'Runx2_area')
        % continue
    end
    file_end_name = 'maxproj.tif';
    proj = imread([paths.inFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(chToUse) file_end_name]);
    proj_blastema = proj(:, myScale.ampPointreal:myScale.endPointrot);
    threshold = graythresh(proj_blastema);
    BW = imbinarize(proj_blastema, threshold);
    area = 0;
    total = 0;
    for row = 1:height(BW)
        for col = 1: width(BW)
            pixel = BW(row, col);
            if pixel ~= 0
                area = area + 1;
                total = total + pixel;
            end
        end
    end
    myScale.Runx2_area = area;
    myScale.Runx2_total_fluorescence = total;
    save(wpathmat, 'myScale');
    close;
end
%%
%%
%%
%% TGMM segmentation starts here
disp(["----------Start TGMM segmentation----------"]);

%% prepare TGMM
% tau_range = [5,6,7,8];
% bckgrd_range = [10,12,15,20];
% tau_range = [5,6,7,8,9,10];
% bckgrd_range = [6,10,12,15,20];
paths.nameFolder=dataFolder_split{end-1};
paths.masterFolder=dataFolder;

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
    datasource   = [paths.masterFolder  'divided_eq' filesep]; % I use equalized images
%     st_dir=dir(strcat(datasource,'fish*','_ch2_hypo.tif')); 
    st_dir=dir(strcat(datasource,'fish',  num2str(fish), '*ray', num2str(scale), '*', num2str(hpp), 'hpa*',...
        '_ch', num2str(ch_nuc), '_hypo.tif')); 
else
    tgmmFolder = 'TGMM_equalized_ch2';
    datasource   = [paths.masterFolder  'equalized' filesep]; % I use equalized images
%     st_dir=dir(strcat(datasource,'fish*','_ch2_hypo.tif')); 
    st_dir=dir(strcat(datasource,'fish',  num2str(fish), '*ray', num2str(scale), '*', num2str(hpp), 'hpa*',...
        '_ch', num2str(ch_nuc), '.tif')); %changed from hypo to none July6
end


tpl=           [TGMMtemplateFolder 'segmentation_package' filesep 'TGMM_configFile.txt'];  
bat_tpl=           [TGMMtemplateFolder 'segmentation_package' filesep 'batch_TGMM.bat'];  

configfolder=  [paths.masterFolder tgmmFolder filesep 'TGMMconfig' filesep];
    


datapathmac  = [paths.masterFolder  tgmmFolder filesep 'data' filesep];
respathmac   = [paths.masterFolder  tgmmFolder filesep 'results' filesep];
    
datapathwin= [winFolder filesep paths.nameFolder filesep tgmmFolder filesep 'data' filesep];
respathwin=  [winFolder filesep paths.nameFolder filesep tgmmFolder filesep 'results' filesep];
    
    
mkdir(configfolder);
mkdir(respathmac);
mkdir(datapathmac);
  

for i=1:numel(st_dir)
    disp(st_dir(i).name);
    prepareTGMM(st_dir(i).name,tpl,configfolder,datasource,datapathmac,respathmac,...
        datapathwin,respathwin,bat_tpl,[paths.masterFolder,tgmmFolder],TGMMbuildFolder...
        ) % background, tau
end

% parameter scanning
% for i=1:numel(st_dir)
%     for tau = tau_range
%         for bckgrd = bckgrd_range
%     disp(st_dir(i).name);
%     prepareTGMM(st_dir(i).name,tpl,configfolder,datasource,datapathmac,respathmac,...
%         datapathwin,respathwin,bat_tpl,[paths.masterFolder,tgmmFolder],TGMMbuildFolder,...
%         bckgrd,tau) % background, tau
%         end
%     end
% end

fclose all;
%% Run TGMM commandline 
paths.masterFolder=dataFolder;

if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end

command = [paths.masterFolder(1:2) ' & cd ' paths.masterFolder tgmmFolder ' & batch_TGMM.bat']; % & here maybe specific for windows
[status,cmdout] = system(command);



%% Parsing TGMM output
%% Organize TGMM output in Matlab structs
% When data sample is tested, this section is skipped as we already provide
% TGMM output in Matlab form.  Users can go to "START ERK  QUANTIFICATION"
% 
paths.masterFolder=dataFolder;

if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end

addpath(genpath(readTGMMFolder))

respath=[dataFolder filesep tgmmFolder filesep 'results' filesep];
objpath=[dataFolder filesep tgmmFolder filesep 'objects' filesep];

mkdir(objpath);

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([respath 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*']);

for i=1:numel(st_dir) 
    disp('-------------------')
    if exist([objpath st_dir(i).name '.mat'],'file')
        % continue
    end
    disp(['Parse tgmm ' st_dir(i).name])
    gme=dir([respath st_dir(i).name filesep 'G*']);
    [~,sort_idx] = sort_nat({gme.date});
    gme = gme(sort_idx);

   if(numel(gme)>0)
     f=dir([respath st_dir(i).name filesep gme(end).name  filesep 'XML_finalResult_lht' filesep 'GME*0000.xml']);

     if(numel(f)>0)

        fxml=([respath st_dir(i).name filesep gme(end).name filesep 'XML_finalResult_lht' filesep f(1).name]);
        obj=readXMLmixtureGaussians(fxml);
        [svList,sizeIm]=readListSupervoxelsFromBinaryFile([fxml(1:end-4) '.svb']);
     end
   end
   save([objpath st_dir(i).name],'obj','svList','sizeIm');
   clear('obj','svList','sizeIm');
end

rmpath(genpath(readTGMMFolder))

%% Measure Erk activity in individual cells
if ~isempty(ch_KTR)

disp('----------Starting measure Erk Activity----------')
if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end
paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.plotFolder=    [paths.masterFolder 'nuclei_example/']; % where the roi is stored
mkdir(paths.plotFolder);

if osteoblast
paths.inFolder=     [paths.masterFolder 'divided_refEq' filesep]; % input folder non-equalized Erk sensor signal
paths.refFolder=    [paths.masterFolder 'divided_eq' filesep]; % input folder equalized nuclear signals
paths.KTRsuffix = ['_ch' num2str(ch_KTR) '_hypo']; % suffix Erk sensor signal stacks
paths.H2Asuffix = ['_ch' num2str(ch_nuc) '_hypo']; % suffix nuclear signal stacks
else
paths.inFolder=     [paths.masterFolder 'rotated' filesep]; % input folder non-equalized Erk sensor signal
paths.refFolder=    [paths.masterFolder 'equalized' filesep]; % input folder equalized nuclear signals
paths.KTRsuffix = ['_ch' num2str(ch_KTR)]; % suffix Erk sensor signal stacks
paths.H2Asuffix = ['_ch' num2str(ch_nuc)]; % suffix nuclear signal stacks
end

paths.objFolder=    [paths.masterFolder 'objects' filesep]; % myScale objects folder
paths.tgmmFolder=   [paths.masterFolder tgmmFolder filesep 'objects' filesep]; % TGMM objects folder


opts=[];
opts.infoName=tgmmFolder; % struct where TGMM information must be stored
opts.saveName='ktr'; % name ktr field
opts.verbose = 0; % display steps processing (0: none; 1: some; 2: all)
opts.checksegmentation = true;
opts.ktrQC = true;

% fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
% scale=''; 
% hpp='';
st_dir=dir([paths.tgmmFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

for i=1:numel(st_dir)

    load([paths.objFolder filesep erase(st_dir(i).name,paths.H2Asuffix)]); %load myScale
    if isfield(myScale,'TGMM_equalized_ch2')&&isfield(myScale.TGMM_equalized_ch2,'ktr')
        % continue;
    end
   disp('------------');
   display([st_dir(i).name]); 
 
   load([st_dir(i).folder filesep st_dir(i).name]); %load TGMM
   load([paths.objFolder filesep erase(st_dir(i).name,paths.H2Asuffix)]); %load myScale
%    load([paths.objFolder filesep 'fish1_ray2_72hpa.mat']); %modified for scanning paramter space

%     myScale.TGMM_equalized_ch2.amputationPoint =  myScale.TGMM_hypo_eq_ch2.amputationPoint;
%     myScale.TGMM_equalized_ch2.endPoint =  myScale.TGMM_hypo_eq_ch2.endPoint;
%    myScale = rmfield(myScale,"TGMM_hypo_eq_ch2");

   [myScale,obj] = assignKTR(myScale,obj,svList,paths,opts);

   save([paths.tgmmFolder filesep st_dir(i).name],'obj','-append');
   save([paths.objFolder filesep erase(st_dir(i).name,paths.H2Asuffix)],'myScale'); 
end

end

%% Measure GEM activity in individual cells
% ch_signal = ch_of_interest(ch_of_interest ~= ch_ref);
if ~isempty(ch_GEM)

disp('----------Starting measure GEM Activity----------')
if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end
paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.plotFolder=    [paths.masterFolder 'nuclei_example/']; % where the roi is stored
mkdir(paths.plotFolder);

if osteoblast
paths.inFolder=     [paths.masterFolder 'divided_refEq/']; % input folder non-equalized Erk sensor signal
paths.refFolder=    [paths.masterFolder 'divided_eq/']; % input folder equalized nuclear signals
paths.GEMsuffix = ['_ch' num2str(ch_GEM) '_hypo']; % suffix Erk sensor signal stacks
paths.H2Asuffix = ['_ch' num2str(ch_nuc) '_hypo']; % suffix nuclear signal stacks
else
paths.inFolder=     [paths.masterFolder 'rotated/']; % input folder non-equalized Erk sensor signal
paths.refFolder=    [paths.masterFolder 'equalized/']; % input folder equalized nuclear signals
paths.GEMsuffix = ['_ch' num2str(ch_GEM)]; % suffix Erk sensor signal stacks
paths.H2Asuffix = ['_ch' num2str(ch_nuc)]; % suffix nuclear signal stacks
end

paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
paths.tgmmFolder=   [paths.masterFolder tgmmFolder '/objects/']; % TGMM objects folder


opts=[];
opts.infoName=tgmmFolder; % struct where TGMM information must be stored
opts.saveName='gem'; % name ktr field
opts.verbose = 0; % display steps processing (0: none; 1: some; 2: all)
if isempty(ch_KTR)
opts.checksegmentation = true;
end

% fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
% scale='2'; 
% hpp='336';
st_dir=dir([paths.tgmmFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

for i=1:numel(st_dir)

   load([paths.objFolder filesep erase(st_dir(i).name,paths.H2Asuffix)]); %load myScale
    if isfield(myScale,'TGMM_equalized_ch2')&&isfield(myScale.TGMM_equalized_ch2,'gem')
        % continue;
    end

   disp('------------');
   display([st_dir(i).name]); 
 
   load([st_dir(i).folder filesep st_dir(i).name]); %load TGMM
%    load([paths.objFolder filesep 'fish1_ray2_72hpa.mat']); %modified for scanning paramter space

%     myScale.TGMM_equalized_ch2.amputationPoint =  myScale.TGMM_hypo_eq_ch2.amputationPoint;
%     myScale.TGMM_equalized_ch2.endPoint =  myScale.TGMM_hypo_eq_ch2.endPoint;
%    myScale = rmfield(myScale,"TGMM_hypo_eq_ch2");

    

   [myScale,obj] = assignGEM(myScale,obj,svList,paths,opts);

   save([paths.tgmmFolder filesep st_dir(i).name],'obj','-append');
   save([paths.objFolder filesep erase(st_dir(i).name,paths.H2Asuffix)],'myScale'); 
end

end

%%
disp('----------segmentation pipeline finished----------')
% toc(tStart)


%% Rotate in xy plane and find amp/end
disp(["----------Start rotation in xy plane and find amp/end points----------"]);
rot_or_not = rot_or_not_post_tgmm;

if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end
paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
paths.roiFolder=    [paths.masterFolder 'rotated/']; % where the roi is stored
paths.finalmaskFolder=    [paths.masterFolder 'finalmask/']; % where the roi is stored

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

% writing options
options=[];
options.overwrite='true';
options.compress='lzw';
options.message = false;

for i=1:numel(st_dir)
%%
    disp('------------');
    display([st_dir(i).name]); 
    
    % load 
    if isfile([paths.roiFolder filesep erase(st_dir(i).name,'.mat') '_ROItight.tif'])
        roipath = [paths.roiFolder filesep erase(st_dir(i).name,'.mat') '_ROItight.tif'];
    elseif isfile([paths.roiFolder filesep erase(st_dir(i).name,'.mat') '_ROI.tif'])
        roipath = [paths.roiFolder filesep erase(st_dir(i).name,'.mat') '_ROI.tif'];
    end
    roiStack = loadtiff([roipath]);
    roiProj = logical(max(roiStack,[],3));

    load([paths.objFolder filesep st_dir(i).name]); %load myScale
    
    cc = myScale.(tgmmFolder).centers;


    % dont rotate if ray is too short % edited 122823
    [coeff,~,latent] = pca(cc(:,[1,2]));
    aspectRatio = sqrt(latent(1)./latent(2));

    if aspectRatio<3 || ~rot_or_not
        th = 0;
        ccrot = cc;
        roiProjRot = roiProj;
    else
        if coeff(2,1)./coeff(1,1)<1 % edited 080124 deal with vertically stretched shape
            f1 = fit(cc(:,1),cc(:,2),'poly1');
            % calculate the slope of this line as a rotation angle
            th = atan(f1.p1);
        else 
            slope_th = coeff(2,2)./coeff(1,2);
            th = atan(slope_th);
        end
        R = rotz(-th*180/pi);
        %     roiProjRot = imrotate(roiProj,th*180/pi);
        [roiProjRot,Rout] = imrotate_tform(roiProj,th*180/pi); % edited 122823
        ccrot = (R*cc')'-[Rout.XWorldLimits(1) Rout.YWorldLimits(1) 0];
    end

%     % dont rotate if ray is too short % edited 121223
%     raylength = max(cc(:,1))-min(cc(:,1));
%     raywidth = max(cc(:,2))-min(cc(:,2));
%     if raylength*0.7<raywidth
%         th = 0;
%     end

    % we rotate the points, see documentation of rotz and Wikipedia
    % "Rotation matrix"
    

    % save finalmask
    saveastiff(im2uint8(roiProjRot),[paths.finalmaskFolder filesep st_dir(i).name(1:end-4) '_ch2maxproj.tif'],options);
    
    % calculate amputation and end points
    roi_ind = find(logical(roiProjRot));
    [roiIdx1,roiIdx2] = ind2sub(size(roiProjRot),roi_ind);

    ampPoint = min(roiIdx2);
    endPoint = max(roiIdx2);

    % update myScale
    myScale.(tgmmFolder).ccrot = ccrot;
    myScale.ampPoint = ampPoint;
    myScale.endPoint = endPoint;
    myScale.pxsize = myScale.metadata{1}.voxel(1);
    myScale.xyRotAngle = -th*180/pi;

    % save myScale
    save([paths.objFolder filesep st_dir(i).name],'myScale'); 

end


%% generate tissue mask from tgmm (xy rot or not)
rot_or_not = rot_or_not_post_tgmm;

if geometry == 1

    disp('----------Calculate geometry masks----------')


if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end

paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
paths.tgmmFolder=   [paths.masterFolder tgmmFolder '/objects/']; % TGMM objects folder
paths.outFolder = [paths.masterFolder filesep '\tissue\'];

if osteoblast
paths.stackFolder =     [paths.masterFolder 'divided_refEq/']; % input folder non-equalized Erk sensor signal
paths.suffix = ['_ch' num2str(ch_KTR) '_hypo']; % suffix Erk sensor signal stacks
paths.H2Asuffix = ['_ch' num2str(ch_nuc) '_hypo']; % suffix nuclear signal stacks
else
paths.stackFolder =     [paths.masterFolder 'rotated/']; % input folder non-equalized Erk sensor signal
paths.suffix = ['_ch' num2str(ch_KTR)]; % suffix Erk sensor signal stacks
paths.H2Asuffix = ['_ch' num2str(ch_nuc)]; % suffix nuclear signal stacks
end


mkdir(paths.outFolder)


% fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
% scale=''; 
% hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

options=[];
options.overwrite='true';
options.compress='lzw';
options.message = false;

for i=1:numel(st_dir)
%
% i = 15;

    if exist([paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_nuclei_eroded.tif'],'file')
        continue
    end

   disp('------------');
   display([st_dir(i).folder filesep st_dir(i).name]);
   load([st_dir(i).folder filesep st_dir(i).name]); %myScale
   load([paths.tgmmFolder filesep st_dir(i).name(1:end-4) paths.H2Asuffix '.mat']); %load TGMM



   
    % generate full mask for entire image
    full_mask = zeros(sizeIm,'logical');
    svList_idx_all = cat(2,obj.svIdx)'+1;
    cell_indices_all = cat(1,svList{svList_idx_all})+1;
    full_mask(cell_indices_all) = 1;
    full_mask = permute(full_mask,[2,1,3]);
    saveastiff(im2uint8(full_mask),[paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_nuclei.tif'],options);


    %% nuclei mask but most outer bound eroded for each nuclei
    full_mask_erode = zeros(sizeIm,'logical');
    for n = 1:numel(obj)
        maskHere = zeros(sizeIm,'logical');
        indicesHere = cat(1,svList{obj(n).svIdx+1})+1;
        maskHere(indicesHere) = 1;
        maskHere_eroded = imerode(maskHere,strel("cube",3));
        % maskHere_eroded = imerode(maskHere,strel("sphere",1));
        full_mask_erode = full_mask_erode|maskHere_eroded;
    end
    full_mask_erode = permute(full_mask_erode,[2,1,3]);
    saveastiff(im2uint8(full_mask_erode),[paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_nuclei_eroded.tif'],options);
    % saveastiff(im2uint8(full_mask_erode),[paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_nuclei_eroded_sphere.tif'],options);

    % labeled image of the nuclei segmentation
    full_mask_labeled = zeros(sizeIm,'double');
    nuclei_labels = cell2mat(...
        arrayfun(@(s,i)i.*ones(size(cat(1,svList{s.svIdx+1}))),...
        obj,(1:numel(obj))',uni=0)...
        );
    full_mask_labeled(cell_indices_all) = nuclei_labels;
    full_mask_labeled = permute(full_mask_labeled,[2,1,3]);
    saveastiff(uint16(full_mask_labeled),[paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_nuclei_labeled.tif'],options);









    %% rotate for tissue mask
    % tic;
    if rot_or_not
        full_mask_rot = [];
        for z = 1:size(full_mask,3)
            imHere = imrotate(full_mask(:,:,z), -myScale.xyRotAngle);
            full_mask_rot = cat(3,full_mask_rot, imHere);
        end
    else
        full_mask_rot = full_mask;
    end
    % toc


    % generate tissue mask
    % tic;
    % tissue_mask = imclose(full_mask_rot,strel('disk',30));
    % clear full_mask_rot
    % tissue_mask = permute(tissue_mask,[1,3,2]);
    % tissue_mask = imclose(tissue_mask,strel('disk',30));
    % tissue_mask = permute(tissue_mask,[1,3,2]);
    % toc
    % tissue_mask_3d = tissue_mask;

    % 2d way of closing
    tic;
    tissue_mask = false(size(full_mask_rot));
    for sl = 1:size(tissue_mask,3) % sl: slice
        tissue_mask(:,:,sl) = imclose(full_mask_rot(:,:,sl),strel('disk',30));
    end
    clear full_mask_rot
    tissue_mask = permute(tissue_mask,[1,3,2]);
    for sl = 1:size(tissue_mask,3) % sl: slice
        tissue_mask(:,:,sl) = imclose(tissue_mask(:,:,sl),strel('disk',30));
    end
    tissue_mask = permute(tissue_mask,[1,3,2]);
    toc

    % figure; sliceViewer(tissue_mask_3d);
    % figure; sliceViewer(tissue_mask);

    tissue_mask = imfill(tissue_mask,"holes");

    saveastiff(im2uint8(tissue_mask),[paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_tgmm.tif'],options);

    % generate z map for tissue mask
    tissue_mask_z_thickness = sum(tissue_mask,3);
    saveastiff(tissue_mask_z_thickness,[paths.outFolder filesep st_dir(i).name(1:end-4) '_ROI_tgmm_z_thickness.tif'],options);

    % calculate crossArea
    crossArea = sum(tissue_mask,[1,3]);

    myScale.crossArea = crossArea;
    % save myScale
    save([paths.objFolder filesep st_dir(i).name],'myScale'); 

    clear full_mask full_mask_erode tissue_mask

end


end


%% Quality control
disp('----------Start QC----------')
if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end
paths=[];
paths.masterFolder=dataFolder; %folder where data is stored
paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
paths.tgmmFolder=   [paths.masterFolder tgmmFolder '/objects/']; % TGMM objects folder
paths.plotFolder=    [paths.masterFolder 'QC/']; % where the roi is stored
mkdir(paths.plotFolder);

% fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
% scale=''; 
% hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);



mkdir([paths.plotFolder filesep 'H2A_by_volume'])
mkdir([paths.plotFolder filesep 'ERKratio_by_ktrAve'])
mkdir([paths.plotFolder filesep 'ktrAve_hist'])
mkdir([paths.plotFolder filesep 'ERKratio_hist'])

for i=1:numel(st_dir)

    disp('------------');
    display([st_dir(i).name]); 
    
    % H2A by volume    
    load([paths.objFolder filesep st_dir(i).name],'myScale'); 

    f = figure('Visible','off');
    x = myScale.(tgmmFolder).volume;
    y = myScale.(tgmmFolder).averageH2ANuc;
    xmax = quantile(x,0.9).*1.2;
    ymax = quantile(y,0.9).*1.2;
%     y(isinf(y)) = 0;
    xplot = x./xmax; yplot = y./ymax;
    scatplot(xplot,yplot,...
        'circles',0.05,100,5,2,4);
    xline(quantile(xplot,0.1),'r')
    xline(quantile(xplot,0.25),'r')
    xline(quantile(xplot,0.9),'r')
    yline(quantile(yplot,0.1),'r')
    yline(quantile(yplot,0.25),'r')
    yline(quantile(yplot,0.9),'r')

    %     xlim([0,4000]);
%     ylim([0,100]);
    xlim([0,1]);
    ylim([0,1]);
    title(st_dir(i).name(1:end-4),'Interpreter','none')
    xlabel('volume')
    ylabel('H2B signal')
    c = colorbar;
    config_plot(f,c);
    xticks([0,0.5,1]);
    xticklabels(arrayfun(@num2str,round(linspace(0,xmax,3),0),'uni',0));
    yticks([0,0.5,1]);
    yticklabels(arrayfun(@num2str,round(linspace(0,ymax,3),0),'uni',0));

    % save
    saveas(f,[paths.plotFolder,filesep,'H2A_by_volume',filesep, st_dir(i).name(1:end-4),'.png']);
    saveas(f,[paths.plotFolder,filesep,'H2A_by_volume',filesep, st_dir(i).name(1:end-4),'.fig']);
    
    close all;

    if ~isempty(ch_KTR)

    % ERK ratio by ktrAve
    ktrNuc = myScale.(tgmmFolder).ktrNuc;
    ktrCyt = myScale.(tgmmFolder).ktrCyt;
    areaNuc = myScale.(tgmmFolder).areaNucleus;
    areaCyt = myScale.(tgmmFolder).areaCyt;
    
    ktrTotal = ktrNuc.*areaNuc+ktrCyt.*areaCyt;
    ktrAve = ktrTotal./(areaNuc+areaCyt);

    myScale.(tgmmFolder).ktrTotal = ktrTotal;
    myScale.(tgmmFolder).ktrAve = ktrAve;

    f = figure('Visible','off');
    xplot = ktrAve;
    yplot = myScale.(tgmmFolder).ktr;
    scatter(xplot,yplot,'.')
    xlim([0,1.2.*quantile(xplot,0.9)])
    ylim([0,2])
    xlabel("average cell KTR")
    ylabel("ERK ratio")
    xline(quantile(xplot,0.25),'r')

    % save
    saveas(f,[paths.plotFolder,filesep,'ERKratio_by_ktrAve',filesep, st_dir(i).name(1:end-4),'.png']);
    saveas(f,[paths.plotFolder,filesep,'ERKratio_by_ktrAve',filesep, st_dir(i).name(1:end-4),'.fig']);
    
    close all;


    % save myScale
    save([paths.objFolder filesep st_dir(i).name],'myScale'); 

    % ktrAve histogram
    f = figure('Visible','off');
    histogram(xplot,'BinLimits',[0,1.2.*quantile(xplot,0.9)]);
    xlabel("average cell KTR")
    xlim([0,inf])

    % save
    saveas(f,[paths.plotFolder,filesep,'ktrAve_hist',filesep, st_dir(i).name(1:end-4),'.png']);
    saveas(f,[paths.plotFolder,filesep,'ktrAve_hist',filesep, st_dir(i).name(1:end-4),'.fig']);
    
    close all;

    % ERK ratio histogram
    f = figure('Visible','off');
    histogram(yplot);
    xlabel("ERK ratio")
    xlim([0,2])
    % save
    saveas(f,[paths.plotFolder,filesep,'ERKratio_hist',filesep, st_dir(i).name(1:end-4),'.png']);
    saveas(f,[paths.plotFolder,filesep,'ERKratio_hist',filesep, st_dir(i).name(1:end-4),'.fig']);
    

    end

    close all;

end
%% plot ERK vs x
% if ~isempty(ch_KTR)
% 
% disp('----------Start plotting results ERK----------')
% if osteoblast
%     tgmmFolder = 'TGMM_hypo_eq_ch2';
% else
%     tgmmFolder = 'TGMM_equalized_ch2';
% end
% paths=[];
% paths.masterFolder=dataFolder; %folder where data is stored
% paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
% paths.tgmmFolder=   [paths.masterFolder tgmmFolder '/objects/']; % TGMM objects folder
% paths.plotFolder=    [paths.masterFolder 'results/']; % where the roi is stored
% mkdir(paths.plotFolder);
% mkdir([paths.plotFolder filesep 'ERK_by_x'])
% 
% % fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
% % scale=''; 
% % hpp='';
% st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
% 
% 
% 
% 
% 
% for i=1:numel(st_dir)
% 
%     display([st_dir(i).name]); 
% 
%     load([paths.objFolder filesep st_dir(i).name]); %load myScale
%     if ~(osteoblast)
%         [myScale,QC_logical] = segmentation_qualitycontrol(myScale,volume = 0.1,averageH2ANuc = 0.1);
%     end
%     ktrAve = myScale.(tgmmFolder).ktrAve;
% 
%     ampPoint =myScale.ampPoint;
%     endPoint =myScale.endPoint;
% 
%     xplot = myScale.(tgmmFolder).ccrot(:,1);
%     yplot = myScale.(tgmmFolder).ktr;
% 
%     if osteoblast
%         trim_logical = xplot>=ampPoint;
%     else
%         trim_logical = xplot>=ampPoint & ktrAve>=quantile(ktrAve,0.25);
%     end
% 
%     trim_logical = true(size(xplot));
% 
% 
%     %plot
%     xplot = xplot(trim_logical)-ampPoint;
%     yplot = yplot(trim_logical);
% 
% 
%     % plot
%     f = figure('Visible','off');
%     cm = colormap(lines(10));
%     scatter(xplot,yplot,'.')
%     ylim([0.4,1.6])
% %     xlim([400,inf])
% 
%     removed = remove_nan_inf([xplot,yplot]);
%     xplot = removed(:,1);
%     yplot = removed(:,2);
%     % fit 
%     fit_results = fit(xplot,yplot,fittype({'x','1'}));
%     hold on;
%     xfit = linspace(min(xplot),max(xplot),100);
%     plot(xfit,fit_results(xfit),linewidth = 2);
%     hold off;
% 
%     % binned average
%     [binned.y,binned.x,binned.sem] = bin_average(xplot,yplot,20);
%     hold on;
%     errorbar(binned.x,binned.y,binned.sem,'-','color',cm(3,:),'LineWidth',2)
%     hold off;
% 
%     xlabel('x')
%     ylabel('ERK')
% 
% %     set(gca, 'XDir','reverse')
% 
%     config_plot(f);
%     exportgraphics(f,[paths.plotFolder filesep 'ERK_by_x' filesep st_dir(i).name(1:end-4) '.png'])
%     close all;
% 
% end
% 
% end
% 
% %% plot GEM vs x
% % ch_signal = ch_of_interest(ch_of_interest ~= ch_ref);
% if ~isempty(ch_GEM)
% 
% disp('----------Start plotting results GEM----------')
% if osteoblast
%     tgmmFolder = 'TGMM_hypo_eq_ch2';
% else
%     tgmmFolder = 'TGMM_equalized_ch2';
% end
% paths=[];
% paths.masterFolder=dataFolder; %folder where data is stored
% paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
% paths.tgmmFolder=   [paths.masterFolder tgmmFolder '/objects/']; % TGMM objects folder
% paths.plotFolder=    [paths.masterFolder 'results/']; % where the roi is stored
% mkdir(paths.plotFolder);
% mkdir([paths.plotFolder filesep 'GEM_by_x'])
% 
% % fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
% % scale=''; 
% % hpp='';
% st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
% 
% 
% 
% 
% 
% for i=1:numel(st_dir)
% 
%     display([st_dir(i).name]); 
% 
%     load([paths.objFolder filesep st_dir(i).name]); %load myScale
%     if ~(osteoblast)
%         [myScale,QC_logical] = segmentation_qualitycontrol(myScale,volume = 0.1,averageH2ANuc = 0.1);
%     end
%     ktrAve = myScale.(tgmmFolder).ktrAve;
% 
%     ampPoint =myScale.ampPoint;
%     endPoint =myScale.endPoint;
% 
%     xplot = myScale.(tgmmFolder).ccrot(:,1);
%     yplot = myScale.(tgmmFolder).gem;
% 
%     if osteoblast
%         trim_logical = xplot>=ampPoint;
%     else
%         trim_logical = xplot>=ampPoint & ktrAve>=quantile(ktrAve,0.25);
%     end
% 
% 
% 
% 
%     %plot
%     xplot = xplot(trim_logical)-ampPoint;
%     yplot = yplot(trim_logical);
% 
% 
%     % plot
%     f = figure('Visible','off');
%     cm = colormap(lines(10));
%     scatter(xplot,yplot,'.')
%     ylim([0 3])
% %     xlim([400,inf])
% 
%     removed = remove_nan_inf([xplot,yplot]);
%     xplot = removed(:,1);
%     yplot = removed(:,2);
%     % fit 
%     fit_results = fit(xplot,yplot,fittype({'x','1'}));
%     hold on;
%     xfit = linspace(min(xplot),max(xplot),100);
%     plot(xfit,fit_results(xfit),linewidth = 2);
%     hold off;
% 
%     % binned average
%     [binned.y,binned.x,binned.sem] = bin_average(xplot,yplot,20);
%     hold on;
%     errorbar(binned.x,binned.y,binned.sem,'-','color',cm(3,:),'LineWidth',2)
%     hold off;
% 
%     xlabel('x')
%     ylabel('GEM')
% 
% %     set(gca, 'XDir','reverse')
% 
%     config_plot(f);
%     exportgraphics(f,[paths.plotFolder filesep 'GEM_by_x' filesep st_dir(i).name(1:end-4) '.png'])
%     close all;
% 
% end
% 
% end
% 
% %% Generate Erk activity maps
% disp('----------Start plotting ERK activity map----------')
% if osteoblast
%     tgmmFolder = 'TGMM_hypo_eq_ch2';
% else
%     tgmmFolder = 'TGMM_equalized_ch2';
% end
% 
% paths=[];
% paths.masterFolder=dataFolder; %folder where data is stored
% paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
% paths.tgmmFolder=   [paths.masterFolder tgmmFolder '/objects/']; % TGMM objects folder
% paths.outFolder=    [paths.masterFolder filesep 'ERK_activity_map/']; % where the roi is stored
% mkdir(paths.outFolder);
% % paths.stackFolder= [paths.masterFolder 'rotated/']; % input folder Erk sensor images
% % paths.suffix='_ch1';
% 
% if osteoblast
% paths.stackFolder =     [paths.masterFolder 'divided_refEq/']; % input folder non-equalized Erk sensor signal
% paths.suffix = ['_ch' num2str(ch_KTR) '_hypo']; % suffix Erk sensor signal stacks
% paths.H2Asuffix = ['_ch' num2str(ch_nuc) '_hypo']; % suffix nuclear signal stacks
% else
% paths.stackFolder =     [paths.masterFolder 'rotated/']; % input folder non-equalized Erk sensor signal
% paths.suffix = ['_ch' num2str(ch_KTR)]; % suffix Erk sensor signal stacks
% paths.H2Asuffix = ['_ch' num2str(ch_nuc)]; % suffix nuclear signal stacks
% end
% 
% 
% % fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
% % scale=''; 
% % hpp='';
% st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);
% 
% opts=[];
% opts.infoName= tgmmFolder; % struct field name for TGMM information
% % opts.statName='stat_hypo_eq_ch2'; % struct field name for Erk information
% opts.adj=[0 254];
% opts.print=true; % print images
% opts.visualType='nuclei'; %  'centers' vs 'nuclei';
% opts.limValue=[0.8 1.2]; % min and max values log(Erk_ratio) for visualization
% opts.imgbox=[4000 4000]; % quantification image size
% opts.dirField='frontDirection'; % struct field name with scale orientation information
% opts.chooseSource=false; % the user can choose manually the source location
% opts.valueToPlot='ktr'; 
% opt.thetaLandmark=0; % rotation to apply
% 
% hppTrues=nan(1,numel(st_dir));
% for i=1:numel(st_dir)
%     load([st_dir(i).folder filesep st_dir(i).name]); %myScale
%     hppTrues(i)=myScale.hppTrue;
% end
% 
% [~,idxHpp]=sort(hppTrues);
% 
% for i=idxHpp
%    obj=[];
%    disp('------------');
%    display([st_dir(i).folder filesep st_dir(i).name]);
%    load([st_dir(i).folder filesep st_dir(i).name]); %myScale
%    load([paths.tgmmFolder filesep st_dir(i).name(1:end-4) paths.H2Asuffix '.mat']); %load TGMM
% 
%    opts.obj=obj;
%    opts.svList=svList;
% 
%     if ~(osteoblast)
%         [myScale,QC_logical] = segmentation_qualitycontrol(myScale);
%         opts.obj = opts.obj(QC_logical);
%     end
% 
%    [myScale] = visualizeERKKTR(myScale,paths,opts);
% %    save([st_dir(i).folder filesep st_dir(i).name],'myScale');
%    close all;
% 
% end
%%


%%
disp('----------segmentation pipeline finished----------')
toc(tStart)
