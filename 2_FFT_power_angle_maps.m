%% These scripts accompany the manuscript:
%  Gulbinaite et al. (2024) "Spatiotemporal resonance in mouse visual
%  cortex" Curr Biol
%
%  MATLAB version used 2019b

%  Requirements (optional): Brewer colormap (https://nl.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps) 
%% ========================================================================
%               Code for FFT ANALYSIS of the pre-processed data
%==========================================================================
% INPUT:
%   - datafiles: *.mat (single-trial data 128px x 128px x time)
% OUTPUT: 
%   - trial-average power and phase maps at the flicker frequency;
%   - single-trial power and phase maps at the flicker frequency (used for
%   computing statistical significance of ssVER maps in 3_Statistics)  

clear all
close all

%% Directories 
rootdir = 'C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Data\'; 
load('C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Code\mymap_30deg.mat')

subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

%% Flicker frequencies used 
freqlist = [2:1:15 18:3:48 52:4:56 64:4:72];

%% Allen map info
Bregref = [10 30];
dx = mymap.bregma(1)-Bregref(1);
dy = mymap.bregma(2)-Bregref(2);

allenmask = imread('C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Code\Allen_CCF_mask.tif','tiff');
allenmask = double(allenmask);

%% 3x3 pixels ROI 
s = size(zeros(128,128));
N = length(s);

% Find the neighbors
[c1{1:N}] = ndgrid(1:3);
c2(1:N) = {2};
offsets = sub2ind(s,c1{:}) - sub2ind(s,c2{:});

%% LOOP over subjects

for subji = 1:length(subjects)
    
    disp(subji)
    cd([rootdir subjects(subji).name '\Ready4preprocessing\'])
    
    % load non-brain pixel mask
    maskfile = dir('Mask*');
    mask0 = imread(maskfile.name,'tiff');
    mask0 = double(mask0);
     
    homedir = [rootdir subjects(subji).name '\Preprocessed_singleTrialGSR_RR_Gaussian\'];
    writedir = [rootdir subjects(subji).name '\Figures\'];
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end
    
    cd(homedir)

    %% Colormap2use
    cmap = colormap(brewermap(256, 'PuOr'));
    cmap = flipud(cmap);
    close all
    %% LOOP over frequencies
    for freqi = 5%1:length(freqlist)
        
        %% Select files with appropriate flicker freq
        filelist = dir(['*' num2str(freqlist(freqi),'%02.f') '_Hz*']);
        filelist = {filelist.name};
        
        %% Load data
       
        for filei = 1:length(filelist)
            disp(filelist{filei})
            filedata = load([ homedir filelist{filei}]);
  
            if filei == 1
                data2analyze = zeros(length(filelist),128*128,length(filedata.time4erp));
            end
            
            % Reapply the mask, because Gaussian filter "extends" the masked area 
            filedata.trialdata_3D = bsxfun(@times, filedata.trialdata_3D, mask0);
            filedata.trialdata_3D = bsxfun(@times, filedata.trialdata_3D, allenmask);
   
            % Reshape into 2D
            data2analyze(filei,:,:) = reshape(filedata.trialdata_3D,128*128,[]);
        end

        %% FFT settings
        num_px = 128*128;
        camera_srate = 150.64;  % empirically defined (see 1_Preprocessing_singleTrial.m file)

        time4ssvep = filedata.time4erp;
        time4fft  = [1000 10000]; % skip the first 1000 ms due to transient response
        tidx      = dsearchn(time4ssvep',time4fft');
        
        resolution = 0.1;
        nFFT = ceil( camera_srate/resolution );
        hz   = camera_srate/2*linspace(0,1,floor(nFFT/2+1));

        skipbins = round(0.5/diff(hz(1:2)));            % skip 0.5 Hz
        numbins  = round(1/diff(hz(1:2))) + skipbins;   % +-1 Hz = 10;
        
        % Initialize matrices to store the results
        fftpower_singletrial_SNR = zeros(num_px, length(hz));
        if freqi == 1
            
            pow_maps_SNR = zeros(length(freqlist),length(freqlist),128*128);
            pow_maps  = zeros(length(freqlist),length(freqlist),128*128);
            angle_maps_av = zeros(length(freqlist),128*128);
            resonance_SNR = zeros(1,length(freqlist));
            resonance_ampl = zeros(1,length(freqlist));
            fftspectrum_SNR = zeros(length(freqlist),length(hz));
            fftspectrum_ampl = zeros(length(freqlist),length(hz));
            
            % Later used for testing statistical significance of ssVER resoponses
            fft_alltrials_amp = zeros(size(data2analyze,1),length(freqlist),length(freqlist),128*128);
            fft_alltrials_snr  = zeros(size(data2analyze,1),length(freqlist),length(freqlist),128*128);
            
        end
        
        %==================================================================
        %%                       FFT 
        %==================================================================
        freq2use = freqlist(freqi);
        freq_idx = dsearchn(hz',freq2use);
        
        % Get ANGLE - trial average
        tmp_av = squeeze(mean(data2analyze(:,:,tidx(1):tidx(2))));
        fft_tmp = fft(tmp_av,nFFT,2)/diff(tidx);
        fft_tmp = fft_tmp(:,floor(1:nFFT/2+1));

        fft_tmp_phase= squeeze(angle(fft_tmp));
        angle_maps_av(freqi,:) = squeeze(fft_tmp_phase(:,freq_idx));

        % Get POWER - single trial
        fft_tmp_av = abs(fft(data2analyze(:,:,tidx(1):tidx(2)),nFFT,3)/diff(tidx)).^2;
        fft_tmp_av = squeeze(fft_tmp_av(:,:,floor(1:nFFT/2+1)));
        
        % Store single-trial data (used later for stat. significance of ssVERs)
        for i = 1:length(freqlist)
            freq2use = freqlist(i);
            freqidx = dsearchn(hz',freq2use);
            if subji == 6 && freqi == 18 % the last trial of subj 6 Rec 5 is missing post-stimulus period and thus was not included
                fft_alltrials_amp(1:9,freqi,i,:) = fft_tmp_av(:,:,freqidx);
            else
                fft_alltrials_amp(:,freqi,i,:) = fft_tmp_av(:,:,freqidx);
            end
            
            % compute SNR
            num = fft_tmp_av(:,:,freqidx);
            denom = mean(fft_tmp_av(:,:,[freqidx-numbins:freqidx-skipbins freqidx+skipbins:freqidx+numbins]),3);
            if subji == 6 && freqi == 18 % the last trial of subj 6 Rec 5 is missing post-stimulus period and thus was not included
                fft_alltrials_snr(1:9,freqi,i,:) = num./denom; 
            else
                fft_alltrials_snr(:,freqi,i,:) = num./denom; % single trials
            end
        end
        
        % average over trials and compute SNR
        fft_tmp_av = squeeze(mean(fft_tmp_av,1));
        fftpower_singletrial = fft_tmp_av;
        
        for hzi = numbins + 1:length(hz) - numbins - 1
            numer = fft_tmp_av(:,hzi);
            denom = mean( fft_tmp_av(:,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) ,2);
            fftpower_singletrial_SNR(:,hzi) = numer./denom; % trial-average
        end

        % store power at all flicker frequencies - used later for
        % noise SNR (H0 or SNR expected when flicker of harmonic 
        % frequenies were not presented 
        for i = 1:length(freqlist)
            freq2use = freqlist(i);
            freqidx = dsearchn(hz',freq2use);
            
            num   = fftpower_singletrial(:,freqidx);
            denom = mean(fftpower_singletrial(:,[freqidx-numbins:freqidx-skipbins freqidx+skipbins:freqidx+numbins]),2);
            pow_maps_SNR(freqi,i,:) = num./denom; % single trials
            pow_maps(freqi,i,:) = num;
        end

        %% Find pixel with the max value at flicker frequency
        [~,maxidx_SNR] = max(squeeze(pow_maps_SNR(freqi,freqi,:))); 
        [~,maxidx] = max(squeeze(pow_maps(freqi,freqi,:)));         
        
        %% Store the FFT results
        resonance_SNR(subji,freqi) = squeeze(nanmean(pow_maps_SNR(freqi,freqi,maxidx_SNR + offsets),3));
        resonance_ampl(subji,freqi) = squeeze(nanmean(pow_maps(freqi,freqi,maxidx + offsets),3));

        % Store spectrum from freq-specific pixel
        fftspectrum_SNR(freqi,:) = nanmean(fftpower_singletrial_SNR(maxidx_SNR + offsets,:));
        fftspectrum_ampl(freqi,:) = nanmean(fftpower_singletrial(maxidx+offsets,:)); 
     
        %% Plotting the results

        clf
        figure(1)
        freq2use = freqlist(freqi);
        freq_idx = dsearchn(hz',freq2use);
        
        subplot(231)
        imagesc(reshape(squeeze(pow_maps_SNR(freqi,freqi,:)),128,128,[]))
        hold on
        % Overlay Allen Map
        plot(Bregref(1),Bregref(2),'b*','markersize',5,'markerfacecolor','w')
        for p = 1:length(mymap.edgeOutline)-1
            plot(mymap.edgeOutline{p}(:, 2)-dx, mymap.edgeOutline{p}(:, 1)-dy,'w-');   
        end
        plot(mymap.bregma(1)-dx,mymap.bregma(2)-dy,'w*','markersize',5,'markerfacecolor','w')
        axis square off
        colorbar
        title([num2str(freq2use) ' Hz'])
        colormap(cmap)

        subplot(232)
        plot(hz,mean(fftpower_singletrial_SNR(maxidx_SNR+offsets,:)),'ko-','Linewidth',1,'MarkerFaceColor','w')
        hold on
        plot(freq2use, mean(pow_maps_SNR(freqi,freqi,maxidx_SNR+offsets),3),'m*','Linewidth',2)
        set(gcf,'color','w')
        ylabel('SNR')
        title(['Power spectrum' num2str(freqlist(freqi)) ' Hz (SNR)'])
        box off
        axis square 
        
        subplot(234)
        imagesc(reshape(squeeze(pow_maps(freqi,freqi,:)),128,128,[]));
        hold on
        % Overlay Allen Map
        plot(Bregref(1),Bregref(2),'b*','markersize',5,'markerfacecolor','w')
        for p = 1:length(mymap.edgeOutline)-1
            plot(mymap.edgeOutline{p}(:, 2)-dx, mymap.edgeOutline{p}(:, 1)-dy,'w-');   
        end
        plot(mymap.bregma(1)-dx,mymap.bregma(2)-dy,'w*','markersize',5,'markerfacecolor','w') 
        axis square off
        colorbar
        title([num2str(freq2use) ' Hz'])

        subplot(235)
        plot(hz,mean(fftpower_singletrial(maxidx+offsets,:)),'ko-','Linewidth',1,'MarkerFaceColor','w')
        hold on
        plot(freq2use, mean(pow_maps(freqi,freqi,maxidx+offsets),3),'m*','Linewidth',2)
        set(gcf,'color','w')
        ylabel('Amplitude (a.u.)')
        title(['Power spectrum' num2str(freqlist(freqi)) ' Hz (ampl) '])
        set(gca,'xlim',[1 55])
        axis square
        box off
     
        ax(2)= subplot(233);
        imagesc(reshape(squeeze(angle_maps_av(freqi,:)),128,128,[]));
        hold on
        plot(Bregref(1),Bregref(2),'b*','markersize',5,'markerfacecolor','w')
        for p = 1:length(mymap.edgeOutline)-1
            plot(mymap.edgeOutline{p}(:, 2)-dx, mymap.edgeOutline{p}(:, 1)-dy,'w-');          
        end
        plot(mymap.bregma(1)-dx,mymap.bregma(2)-dy,'w*','markersize',5,'markerfacecolor','w')
        axis square off
        colormap(ax(2),twilight)
        title('Trial-average angle')
        colorbar
        set(gca,'clim',[-pi pi]);
        
        %% save figures
        figure_name = ([ writedir num2str(freqlist(freqi),'%02.f') '_Hz' ]);
        saveas(gcf, figure_name, 'jpeg')

    end
    %% Save FFT results
    % long version
    outfile = [writedir '\FFT_results_peakmaps_withSpectrum.mat'];
    save(outfile,'pow_maps_SNR','pow_maps','fftspectrum_SNR','fftspectrum_ampl','angle_maps_av','fft_alltrials_amp','fft_alltrials_snr','resonance_SNR','resonance_ampl','hz','resolution','skipbins','numbins')

    % short version
    outfile = [writedir '\FFT_results_peakmaps_short_withSpectrum.mat'];
    save(outfile,'pow_maps_SNR','hz','pow_maps','fftspectrum_SNR','fftspectrum_ampl','angle_maps_av','resonance_SNR','resonance_ampl','resolution','skipbins','numbins') 
end

