
%% These scripts accompany the manuscript:
%  Gulbinaite et al. (2024) "Spatiotemporal resonance in mouse visual
%  cortex" Curr Biol

% MATLAB version used 2019b
%% Code for PREPROCESSING the raw data (responses to 2-72 Hz flicker stimulation)

% PART 1: 
% -------------------------------------------------------------------------
%   Requires: 
%   - MATLAB Image Processing Toolbox and Signal Processing Toolbox
%   - Additional files in the folder: abfload.m; findpeaks.m; imreadallraw.m
%
%   INPUT: 
%        - datafiles: *.raw
%        - stimulus onset/offset timestamp files: *.abf
%        - stimulus files (flicker freq order during the experiment): *Rec.mat
%        - transformation matrix used to align across recordings: *adjusted.mat
%        ( done in 0_Coregistration_across_recordings.m )
%        - matrix marking pixels containing brain data (mask): Mask*.mat  

%   OUTPUT: single-trial *.mat files pixels (only-brain pixels) x time 
%           To go back to a full 128px x 128px x time matrix: 
% 
%           trialdata_3D = zeros(128*128, length(time4erp)); 
%           trialdata_3D(NonZeroPixelsIndex,:) = trialdata;
%           trialdata_3D = reshape(trialdata_3D,128,128,[]);

% PART 2: Robust regression - removes slow haemodynamic trend
% -------------------------------------------------------------------------
% Requires: NoiseTools toolbox: http://audition.ens.fr/adc/NoiseTools/ 


% PART 3: Spatial Gaussian smoothing
% -------------------------------------------------------------------------

% =========================================================================
%%       PART 1:  Load, align between recordings, GSR, baseline-correct
%==========================================================================
%%

clear all
close all

rootdir = '...\Spatiotemporal_resonance_DATA\Data\Flicker'; % add path to the data
subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..','Focused_HD_videos'})) = [];  %remove . and ..

%% LOOP over subjects
for subji = 1:length(subjects)
    disp(subji)
    tic
    homedir = [rootdir subjects(subji).name '\Ready4preprocessing\'];
    writedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial\'];
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end

    %% Find all the necessary filetypes
    cd(homedir)
    datafiles        = dir('*.raw');
    markerfiles      = dir('*.abf');
    matfiles         = dir('*Rec.mat');
    adjustment_files = dir('*adjusted.mat');
    adjustment_files = {adjustment_files.name};
    
    maskfile = dir('Mask_*');
    %% LOOP over files
    for filei = 1:length(markerfiles)
        
        clear stimMarks
        
        %% Load stimulus markers and find indices based on camera time
        % marker_data(:,1) - stimulus waveform
        % marker_data(:,2) - behavioral camera markers
        % marker_data(:,3) - camera clock
        % marker_data(:,4) - camera start/stop
        
        marker_data = abfload([homedir markerfiles(filei).name],'start',1,'stop','e');
        
        num_fr = str2num(markerfiles(filei).name(22:26));
        num_px = 128*128;
        marker_srate = 2000; % 2 kHz
        time = 0:1/marker_srate:(length(marker_data(:,1))/marker_srate)-1/marker_srate;

        ntrials = 30;
        
        % Find stimulus peaks
        %-----------------------------------------------------------------
        [pks, locs] = findpeaks(squeeze(marker_data(:,1)),'MinPeakProminence',0.02);
        
        % Remove the second of the two nearby peaks
        badies = find(locs(2:end)-locs(1:end-1)<10000)+1;
        
        locs(badies) = [];
        pks(badies)  = [];
        
        % plot2check
%         figure
%         plot(time,marker_data(:,1),'b-','Linewidth',1)
%         hold on
%         for i = 1:length(pks)
%             plot(time(locs(i)),pks(i),'r*','Linewidth',1,'Markersize',8)
%         end

        % find camera times for the beginning and end of the recording
        %-----------------------------------------------------------------
        tmpidx = find(marker_data(:,4)>3);
        cameraStartStop = [tmpidx(1) tmpidx(end)];
        
        % plot2check
%         figure
%         plot(marker_data(:,4),'bo-','Linewidth',1)
%         hold on
%         plot([cameraStartStop(1) cameraStartStop(1)],[0 4],'r-')
%         plot([cameraStartStop(2) cameraStartStop(2)],[0 4],'r-')

        % find camera clock indices in abf file that correspond to camera frames
        %---------------------------------------------------------------------
        tmpidx = find(marker_data(:,3)>0.5);
        badies = find(tmpidx(2:end)-tmpidx(1:end-1)<2); % remove if two consequtive points
        tmpidx(badies)=[];
        tmpidx(tmpidx<cameraStartStop(1) | tmpidx>cameraStartStop(2))= [];
        
        camera_idx = tmpidx(1:num_fr); % camera frames as indices in abf file
        
        % find camera frames that denote stimulus onset
        stim_idx = dsearchn(camera_idx, locs);
        
        % Display error message if not all trial onsets were found
        if ntrials ~=length(stim_idx)
            disp('Error: Not all trials were found')
        else
            disp('Success!')
        end
        
        %% Calculate empirical camera sampling rate
        
        % empirical camera sampling rate
%         camera_srate = 1./(mean(diff(camera_idx))/marker_srate);
        
        % Number of frames between each camera snapshot
        % figure
        % hist(diff(camera_idx))
        
        %% Load mask matrix (1 = pixels with the brain data)
        % Same mask for all the recordings from that animal
       
        disp(maskfile)
        mask0 = imread(maskfile.name,'tiff');
        NonZeroPixelsIndex = find(reshape(mask0,128*128,1) == 1);
 
        %% Load the data
        
        data_raw = imreadallraw([homedir markerfiles(filei).name(1:end-4) '.raw'],128,128,num_fr,'*uint16');
        
        % Select transformation matrix file for alignement between recordings
        if filei == 1 || filei == 2
            load([homedir adjustment_files{1}])
            disp(adjustment_files{1})
        elseif filei == 3 || filei == 4
            clear mytform
            load([homedir adjustment_files{2}])
            disp(adjustment_files{2})
        elseif filei == 5 || filei == 6
            clear mytform
            load([homedir adjustment_files{3}])
            disp(adjustment_files{3})
        elseif filei == 7 || filei == 8
            clear mytform
            load([homedir adjustment_files{4}])
            disp(adjustment_files{4})
        elseif filei == 9 || filei == 10
            clear mytform
            load([homedir adjustment_files{5}])
            disp(adjustment_files{5})
        elseif filei == 11 || filei == 12
            clear mytform
            load([homedir adjustment_files{5}])
            disp(adjustment_files{5})
        end

        % Apply transformation matrix
        data_all = zeros(128,128,size(data_raw,3));
        img = zeros(128,128);
        
        parfor ii = 1:size(data_raw,3) 
            data_all(:,:,ii) = imwarp(data_raw(:,:,ii),mytform,'outputview',imref2d(size(img)));
        end
       
        clear data_raw
        
        % double precision
        data_all = double(data_all);

        %% Settings 
        
        camera_srate = 150;
        start_epoch = 2; 
        stop_epoch  = 15; 
        prestim     = start_epoch*camera_srate; 
        nsignal     = stop_epoch*camera_srate; 
        time4erp    = -start_epoch*1000:1/camera_srate*1000:stop_epoch*1000; % time vector
        
        % baseline for computing dF/F0
        baseline = [-1000 -100];
        baseline_idx = dsearchn(time4erp',baseline');
        
        %% Load stimulus condition datafile
        load([homedir matfiles(filei).name])
        
        if str2num(matfiles(filei).name(2)) == 1
            stimMarks = frex(1:ntrials);
        elseif str2num(matfiles(filei).name(2)) == 2
            stimMarks = frex(ntrials+1:end);
        end
        
        freqlist = unique(frex);
        
        %% Loop over trials: (1) perform Global Signal Regression (GSR) (2) Baseline-correct the signal
        
        for triali = 1:ntrials
            
            % get single trial around the marker time
            trialdata = data_all(:,:,stim_idx(triali)-prestim:stim_idx(triali) + nsignal);
            
            % keep only NonZero pixels for computational efficiency
            trialdata = reshape(trialdata,num_px, prestim + nsignal + 1);
            trialdata = trialdata(NonZeroPixelsIndex,:);
            f0 = mean((trialdata(:,baseline_idx(1):baseline_idx(2))),2); 

            %% apply Global signal regression (GSR)
            trialdata = GSR_NonZeroPixels(trialdata)+ mean(f0);

            %% Baseline-correct the signal on a single-trial level
            f0 = mean((trialdata(:,baseline_idx(1):baseline_idx(2))),2); 
            trialdata = 100*bsxfun(@rdivide,bsxfun(@minus,squeeze(trialdata),f0),f0);

            %% Save data as single trials
            trialdata = single(trialdata);

            frequency = stimMarks(triali);
            outfilename = [ writedir datafiles(filei).name(1:2) '_' num2str(triali,'%02.f') 'trial_' num2str(stimMarks(triali),'%02.f') '_Hz.mat'];
            
            save(outfilename,'trialdata','frequency','NonZeroPixelsIndex','time4erp','nsignal','prestim','freqlist','camera_srate','maskfile','baseline','-v7.3');

        end % END Loop over trials     
        
    end % END Loop over recordings
    
end % END Loop over subjects

%==========================================================================
%%                  PART 2: Perform Robust Regression
%==========================================================================
% Requires NoiseTools toolbox: http://audition.ens.fr/adc/NoiseTools/ 

% NOTE: Takes a long time to run and can be skipped, because it is used to 
% remove slow (<1 Hz) trends and is a good practice.
% However, all the subsequent analyses are focused on flicker responses and 
% the slowest flicker frequency used in the study was 2 Hz.

clear all
close all

rootdir = '...\Spatiotemporal_resonance_DATA\Data\Flicker\'; % add path to the data

subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..
%% LOOP over subjects
for subji = 1%:length(subjects)
    disp(subji)
    
    homedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial\'];
    writedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial\'];
    cd(homedir)

    %% List of files
    filelist = dir('*Hz.mat');
    filelist = {filelist.name};
    
    %% Robust regression
    parfor triali = 1:length(filelist) % LOOP over trials
        filename = [ homedir filelist{triali}];
        outfilename = [writedir filelist{triali}(1:end-4) '_withRR.mat'];
        
        % Robust regression
        robust_regress_SSVEP(triali, filename, outfilename);
    end
end

% =========================================================================
%%                  PART 3:  Apply spatial Gaussian smoothing
%==========================================================================
clear all
close all

rootdir = '...\Spatiotemporal_resonance_DATA\Data\Flicker\'; % add path to the data

subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..','Focused_HD_videos'})) = [];  %remove . and ..

%% LOOP over subjects
for subji = 1:length(subjects)
    disp(subji)
   
    homedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial\'];
    writedir = [rootdir subjects(subji).name '\Preprocessed_singleTrialGSR_RR_Gaussian\']; 
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end

    cd(homedir)
    filelist = dir('*withRR.mat');
    filelist = {filelist.name};
    
    %% Loop over trials
    
    for triali = 1:length(filelist)
        
        outfilename = [writedir filelist{triali}(1:end-4) '_gauss.mat'];
        if exist(outfilename,'file'); continue; end
        
        clear trialdata
        disp(filelist{triali});
        load(filelist{triali})
        
        trialdata_3D = zeros(128*128, length(time4erp)); % initialize each time to clean up the temp matrix
        trialdata_3D(NonZeroPixelsIndex,:) = trialdata;
        trialdata_3D = reshape(trialdata_3D,128,128,[]);
        %% Apply Gaussian
        hsize = 8;
        trialdata_3D = gaussianfilter(trialdata_3D, hsize, 2);       
        %% save
        trialdata_3D = single(trialdata_3D);
        save(outfilename,'trialdata_3D','NonZeroPixelsIndex','time4erp','hsize')

    end
end

