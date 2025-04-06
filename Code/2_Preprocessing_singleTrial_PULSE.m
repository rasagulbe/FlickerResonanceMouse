
%% These scripts accompany the manuscript:
%  Gulbinaite et al. (2024) "Spatiotemporal resonance in mouse visual
%  cortex" Curr Biol

% MATLAB version used 2019b

%% Code for PREPROCESSING the raw data (responses to single-pulse stimuli)

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

rootdir = '...\Spatiotemporal_resonance_DATA\Data\Pulse\'; % add path to the data
subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..','Focused_HD_videos'})) = [];  %remove . and ..

%% LOOP over subjects
for subji = 1:length(subjects)

    disp(subji)
    tic
    homedir = [rootdir subjects(subji).name '\'];
    writedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial_ERP\'];
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end

    %% Find all the necessary filetypes
    cd(homedir)
   
    datafiles        = dir('*.raw');
    adjustment_files = dir('*adjusted.mat');
    adjustment_files = {adjustment_files.name};
    markerfiles      = dir('*.abf');
    
    maskfile = dir('Mask_*');

    %% LOOP over files
    for filei = 1:length(markerfiles)
        clear stimMarks
        
        %% Load stimulus markers and find indices based on camera time
        % marker_data(:,1) - stimulus
        % marker_data(:,2) - behavioral camera markers
        % marker_data(:,3) - camera clock
        % marker_data(:,4) - camera start/stop
        
        marker_data = abfload([homedir markerfiles(filei).name],'start',1,'stop','e');
        
        num_fr = str2num(markerfiles(filei).name(17:22)); % pulse

        if isempty(num_fr)
            error('Number of frames undefined!')
            break;
        end
        num_px = 128*128;
        marker_srate = 20000; % 2 kHz
        time = 0:1/marker_srate:(length(marker_data(:,1))/marker_srate)-1/marker_srate;
        
        ntrials = 60;
        
        % Find stimulus peaks
        %-----------------------------------------------------------------
        [pks, locs] = findpeaks(squeeze(marker_data(:,1)),'MinPeakProminence',0.015);
        
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
%                 cameraStartStop(1) = 21960000; % for 7th subject, ERP
        
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
%         tmpidx = tmpidx - cameraStartStop(1); % for subj7
        
        camera_idx = tmpidx; % camera frames as indices in abf file
        
        % find camera frames that denote stimulus onset
        stim_idx = dsearchn(camera_idx, locs);       
%         stim_idx(1:3) = []; % for subj7

        % Display error message if not all trial onsets were found
        if ntrials ~=length(stim_idx)
            disp('Error: Not all trials were found')
        else
            disp('Success!')
        end

        %% Load mask matrix (1 = pixels with the brain data)
        % Same mask for all the recordings from that animal

        disp(maskfile)
        mask0 = imread(maskfile.name,'tiff');
        NonZeroPixelsIndex = find(reshape(mask0,128*128,1) == 1);
 
        %% Load the data
        data_raw = imreadallraw([homedir markerfiles(filei).name(1:end-4) '.raw'],128,128,num_fr,'*uint16');
        
        % Apply transformation for alignement between recordings
        load([homedir adjustment_files{1}])
        disp(adjustment_files{1})

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
        stop_epoch  = 9;
        prestim     = start_epoch*camera_srate;
        nsignal     = stop_epoch*camera_srate;
        time4erp    = -start_epoch*1000:1/camera_srate*1000:stop_epoch*1000; % time vector
        
        % baseline for computing dF/F0
        baseline = [-1000 -100];
        baseline_idx = dsearchn(time4erp',baseline');

        %% Loop over trials: (1) perform Global Signal Regression (GSR) (2) Baseline-correct the signal

        for triali = 1:ntrials
            
            % get single trial around the marker time
            trialdata = data_all(:,:,stim_idx(triali)-prestim:stim_idx(triali) + nsignal);
            
            % keep only NonZero pixels for computational efficiency
            trialdata = reshape(trialdata,num_px, prestim+nsignal+1);
            trialdata = trialdata(NonZeroPixelsIndex,:);
            f0 = mean((trialdata(:,baseline_idx(1):baseline_idx(2))),2);

            %% apply Global signal regression (GSR)
            trialdata = GSR_NonZeroPixels(trialdata)+ mean(f0);
            
            %% Baseline-correct the signal on a single-trial level
            f0 = mean((trialdata(:,baseline_idx(1):baseline_idx(2))),2);
            trialdata = 100*bsxfun(@rdivide,bsxfun(@minus,squeeze(trialdata),f0),f0);
         
            %% Save
            trialdata = single(trialdata);
            outfilename = [ writedir num2str(triali,'%02.f') '_trial'];
            
            save(outfilename,'trialdata','NonZeroPixelsIndex','time4erp','nsignal','prestim','camera_srate','start_epoch','stop_epoch','maskfile','-v7.3');
            
        end % END Loop over trials     
        
    end % END Loop over recordings   
    
end % END Loop over subjects

%==========================================================================
%%                  PART 2: Perform Robust Regression
%==========================================================================
% Requires NoiseTools toolbox: http://audition.ens.fr/adc/NoiseTools/ 

clear all
close all

rootdir = '...\Spatiotemporal_resonance_DATA\Data\Pulse\'; % add path to the data

subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..
%% LOOP over subjects
for subji = 1:length(subjects)

    disp(subji)
    
    homedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial_ERP\'];
    writedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial_ERP\'];
    cd(homedir)

    cd(homedir)
    filelist = dir('*trial.mat');
    filelist = {filelist.name};
    
    %% Settings for Robust Regression
    time2analyze = [-2 4]; 
    camera_srate = 150;
    time4erp    = time2analyze(1)*1000:1/camera_srate*1000:time2analyze(2)*1000; % time vector

    % ERP mask
    start_mask   = round(0*camera_srate);    
    end_mask     = round(0.7*camera_srate);  

    % Create temporal mask for the trial (same for all trials)
    trial_temporal_mask = ones(1,length(time4erp));  % array to store the temporal mask to mask-out all trial onsets
    
    zerodix = dsearchn(time4erp',0);
    mask_startidx  = zerodix; 
    mask_stopidx   = zerodix + end_mask;

    trial_temporal_mask(mask_startidx:mask_stopidx) = 0;
    %% LOOP over files
    tic
    parfor triali = 1:length(filelist)
        
        outfilename = [writedir filelist{triali}(1:end-4) '_withRR.mat'];
        if exist(outfilename,'file')
            continue
        end
        
        filename = [ homedir filelist{triali}];
 
        robust_regress_ERP(triali, filename, outfilename, trial_temporal_mask,time2analyze);
    end % END loop over TRIALS
    toc
end % END loop over SUBJECTS
        
% =========================================================================
%%                  PART 3:  Apply spatial Gaussian smoothing
%==========================================================================

clear all
close all

rootdir = 'D:\Spatiotemporal_resonance_DATA\Data\Pulse\'; % add path to the data

subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..','Focused_HD_videos'})) = [];  %remove . and ..

%% LOOP over subjects
for subji = 1:length(subjects)

    %% Mask file
    cd([rootdir subjects(subji).name])
    maskfile = dir('Mask*');
    mask0 = imread(maskfile.name,'tiff');
    NonZeroPixelsIndex = find(reshape(mask0,128*128,1) == 1);
    %% Datafiles
    
    disp(subji)
   
    homedir = [rootdir subjects(subji).name '\Preprocessed_singleTrial_ERP\'];
    writedir = [rootdir subjects(subji).name '\Preprocessed_singleTrialGSR_RR_Gaussian\']; 
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end

    cd(homedir)
    filelist = dir('*withRR.mat');
    filelist = {filelist.name};
  
    %% Loop over trials
    
    for triali = 1%:length(filelist)
        
        outfilename = [writedir filelist{triali}(1:end-4) '_gauss.mat'];
        if exist(outfilename,'file'); continue; end
        
        clear trialdata
        disp(filelist{triali});
        load(filelist{triali})
        
        % Restore 128 x 128 dimensions
        trialdata_3D = zeros(128*128, length(time4erp)); % initialize it each time to clean up the temp matrix
        trialdata_3D(NonZeroPixelsIndex,:) = trialdata;
        trialdata_3D = reshape(trialdata_3D,128,128,[]);
        %% Apply Gaussian
        hsize = 8;
        trialdata_3D = gaussianfilter(trialdata_3D, hsize, 2);  
        
        %% save
        trialdata_3D = single(trialdata_3D);
        save(outfilename,'trialdata_3D','NonZeroPixelsIndex','time4erp','hsize')
    
    end % END loop over TRIALS
end % END loop over SUBJECTS



