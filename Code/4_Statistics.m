%% These scripts accompany the manuscript:
%  Gulbinaite et al. (2024) "Spatiotemporal resonance in mouse visual
%  cortex" Curr Biol
%
%  MATLAB version used 2019b
%% Code for STATISTICAL significance of flicker responses 

% PART 1: Testing statistically significant responses to flicker in V1 
% -------------------------------------------------------------------------
%  INPUT: Files containing FFT results (from 2_Preprocessing_singleTrial.m)  
%  OUTPUT: t-values

% PART 2: Statistical significance of flicker responses at the cortical map 
% level 
% -------------------------------------------------------------------------
%  INPUT: Files containing FFT results (from 2_Preprocessing_singleTrial.m)  
%  OUTPUT: t-value maps

% NOTES: The code adapted from Pernet et al. (2015) J Neurosci Methods
% "Cluster-based computational methods for mass univariate analyses
% of event-related brain potentials/fields: A simulation study" 

%% ========================================================================
%    PART 1:  Testing statistically significant responses to flicker in V1 
%% ========================================================================
set(0,'DefaultFigureWindowStyle','normal')

clear all
close all
%% Directories 

rootdir = 'C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Data\'; 
subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

%% Flicker frequencies used 

freqlist = [2:1:15 18:3:48 52:4:56 64:4:72];
%% Settings 

s = size(zeros(128,128));
N = length(s);
[c1{1:N}]=ndgrid(1:3);
c2(1:N)={2};
offsets = sub2ind(s,c1{:}) - sub2ind(s,c2{:});

%  Matrix of harmonics
harm_mat = zeros(length(freqlist),30);
for freqi = 1:length(freqlist)
    tmp = freqlist(freqi):freqlist(freqi):80;
    harm_mat(freqi,1:length(tmp)) = tmp;
end

%% LOOP over subjects
for subji = 1:length(subjects)
    
    disp(subji)
    homedir = [rootdir subjects(subji).name '\Figures\']; 

    %% Initialize variables2store
    
    if subji == 1
        percentile_95_SNR = zeros(length(subjects),length(freqlist));
        percentile_95_ampl = zeros(length(subjects),length(freqlist));

        resonance_SNR_all = zeros(length(subjects),length(freqlist));
        resonance_ampl_all = zeros(length(subjects),length(freqlist)); 
    end
    %% Load the data
    filename = [homedir 'FFT_results_peakmaps_withSpectrum.mat'];
    load(filename); 

    %% Store resonance data per subject
    resonance_SNR_all(subji,:) = resonance_SNR(subji,:);
    resonance_ampl_all(subji,:) = resonance_ampl(subji,:);

    %% LOOP over frequencies and get 95% percentile of noise
    alpha = 0.05;
    z_val = norminv(1-alpha./2);
    num_trialsused =[];
    for freqi = 1:length(freqlist)

        % Find peak pixel
        [~,maxidx_SNR] = max(squeeze(pow_maps_SNR(freqi,freqi,:))); 

        % Exclude trials that contain harmonically related stimulus freqs
        freqs2use = freqlist;
        idx2remove = [];
        for i = 1:length(freqlist)
            if ismember(freqlist(freqi),harm_mat(i,:) )
               idx2remove = [idx2remove i];
            end
        end
        freqs2use(idx2remove) = [];
        freqs2use_idx = dsearchn(freqlist',freqs2use');

        % H0: Noise trials (SNR)
        noise_trials_snr = squeeze(mean(fft_alltrials_snr(:,freqs2use_idx,freqi,maxidx_SNR+offsets),4));
        noise_trials_snr = reshape(noise_trials_snr,size(noise_trials_snr,1)*size(noise_trials_snr,2),1);

        % H0: Noise trials (amplitude)
        noise_trials_ampl = squeeze(mean(fft_alltrials_amp(:,freqs2use_idx,freqi,maxidx_SNR+offsets),4));
        noise_trials_ampl = reshape(noise_trials_ampl,size(noise_trials_ampl,1)*size(noise_trials_ampl,2),1);
        
        % 95% percentile
        percentile_95_SNR(subji,freqi) = prctile(noise_trials_snr,100-(100*alpha));
        percentile_95_ampl(subji,freqi) = prctile(noise_trials_ampl,100-(100*alpha));
    end

    %% PLOT
    figure
    plot(freqlist,squeeze(resonance_SNR(subji,:)),'k-','Linewidth',2)
    hold on
    plot(freqlist, percentile_95_SNR(subji,:),'k:','Linewidth',2)
    legend({'ssVER observed' 'H0 hypothesis'})
    legend box off
    set(gca,'xlim',[2 72],'xscale','log', 'xtick',round(logspace(log10(freqlist(1)),log10(freqlist(end)),6)))
    set(gcf,'color','w')
    axis square
    box off
    xlabel('Flicker frequency (Hz)')
    ylabel('ssVER power (SNR)')
    title(['Animal ' num2str(subji)])
 
end % END loop over subjects

%% PLOT subject-average results
ns=7;

figure
shadedErrorBar(freqlist,mean(resonance_SNR_all),std(resonance_SNR_all)/sqrt(ns),'lineProps',{'mo-','MarkerFaceColor','w','MarkerSize',8});
hold on
shadedErrorBar(freqlist,mean(percentile_95_SNR),std(percentile_95_SNR)/sqrt(ns),'lineProps',{'ko-','MarkerFaceColor','w','MarkerSize',8});

set(gca,'xlim',[2 80],'xscale','log', 'xtick',round(logspace(log10(freqlist(1)),log10(freqlist(end)),6)))
legend({'Subject average' 'H0 hypothesis' })
legend boxoff
set(gcf,'color','w')
pbaspect([2 1.5  1])
box off

%% ========================================================================
%    PART 2:  Statistically significant responses to flicker (cortical-map
%    level)
%% ========================================================================

set(0,'DefaultFigureWindowStyle','docked')

clear all
close all

%% Directories 
rootdir = 'C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Data\'; 
subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

%% Load Allen CCF and mask
load('C:\Users\abete\Desktop\Code\mymap_30deg.mat')

allenmask = imread('C:\Users\abete\Desktop\Code\Allen_CCF_mask.tif','tiff');
allenmask = double(allenmask);

% Allen map info
Bregref = [10 30];
dx = mymap.bregma(1)-Bregref(1); % -5-Bregref(1);
dy = mymap.bregma(2)-Bregref(2);

%% Flicker frequencies used 

freqlist = [2:1:15 18:3:48 52:4:56 64:4:72];
%% LOOP over subjects

for subji = 1:length(subjects)
    disp(subji)
    
    homedir = [rootdir subjects(subji).name '\Figures\'];
    writedir = [rootdir subjects(subji).name '\Stats\'];
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end
                                           
    cd(homedir)
    filelist = dir('FFT_results_peakmaps_withSpectrum.mat');

    filelist = {filelist.name};
    load(filelist{1})

    %%  Matrix of harmonics 
    % Is is important to not include trials in which harmonically-related
    % flicker frequencies were used (see the paper)
    
    harm_mat = zeros(length(freqlist),30);
    for freqi = 1: length(freqlist)
        tmp = freqlist(freqi):freqlist(freqi):80;
        harm_mat(freqi,1:length(tmp)) = tmp;
    end
    
    %----------------------------------------------------------------------
    %% LOOP over flicker frequencies and store thresholded tvalue maps 
    %----------------------------------------------------------------------
    [real_t_all,real_t_thresh_bootc_all] = deal(zeros(length(freqlist),128,128));
    
    for freqi = 1:length(freqlist)

        %% z-score single trial maps - SIGNAL
        %==================================================================
        % mean and over pxels - SNR
        mean_pow_real = squeeze(nanmean(fft_alltrials_snr(:,freqi,freqi,:),4));
        std_pow_real  = squeeze(nanstd(fft_alltrials_snr(:,freqi,freqi,:),[],4));
        signal_trials = bsxfun(@rdivide, bsxfun(@minus, squeeze(fft_alltrials_snr(:,freqi,freqi,:)), repmat(mean_pow_real,1,128*128)),repmat(std_pow_real,1,128*128));
    
        %% z-score single trial maps of NOISE
        % fft_alltrials_snr(:,freqi,i,:): 2nd dimension - flicker frequency on those trials
        % 3d dimension - reponse at all the other frequencies; 
        %==================================================================
        if freqi == 1
            mean_pow_real = squeeze(nanmean(fft_alltrials_snr(:,freqi+1,freqi,:),4));
            std_pow_real  = squeeze(nanstd(fft_alltrials_snr(:,freqi+1,freqi,:),[],4));
            noise_trials = bsxfun(@rdivide, bsxfun(@minus, squeeze(fft_alltrials_snr(:,freqi+1,freqi,:)), repmat(mean_pow_real,1,128*128)),repmat(std_pow_real,1,128*128));
            
        elseif freqi == length(freqlist)
            mean_pow_real = squeeze(nanmean(fft_alltrials_snr(:,freqi-1,freqi,:),4));
            std_pow_real  = squeeze(nanstd(fft_alltrials_snr(:,freqi-1,freqi,:),[],4));
            noise_trials = bsxfun(@rdivide, bsxfun(@minus, squeeze(fft_alltrials_snr(:,freqi-1,freqi,:)), repmat(mean_pow_real,1,128*128)),repmat(std_pow_real,1,128*128));
            
        else  
            mean_pow_real = squeeze(nanmean(fft_alltrials_snr(:,freqi-1,freqi,:),4));
            std_pow_real  = squeeze(nanstd(fft_alltrials_snr(:,freqi-1,freqi,:),[],4));
            noise_trials_low = bsxfun(@rdivide, bsxfun(@minus, squeeze(fft_alltrials_snr(:,freqi-1,freqi,:)), repmat(mean_pow_real,1,128*128)),repmat(std_pow_real,1,128*128));
            
            mean_pow_real = squeeze(nanmean(fft_alltrials_snr(:,freqi+1,freqi,:),4));
            std_pow_real  = squeeze(nanstd(fft_alltrials_snr(:,freqi+1,freqi,:),[],4));
            noise_trials_hi = bsxfun(@rdivide, bsxfun(@minus, squeeze(fft_alltrials_snr(:,freqi+1,freqi,:)), repmat(mean_pow_real,1,128*128)),repmat(std_pow_real,1,128*128));
         
            noise_trials = cat(1,noise_trials_low,noise_trials_hi);
            clear noise_trials_low noise_trials_hi
        end
        
        % Account for one missing trial in freqi=18 cond for subj=6
        if subji == 6 && freqi==17
            noise_trials(20,:,:) = [];
          
        elseif subji == 6 && freqi==18
            signal_trials(10,:,:) = [];          
        elseif subji == 6 && freqi==19
           noise_trials(10,:,:) = [];
        end
      
        %% Compute REAL t-value map
        %==================================================================
      
        pval_px = 0.001;              % pixel p-value threshold
        pval_cluster = 0.01;          % cluster correction p-value threshold
        nTrials = size(signal_trials,1);
        zval = norminv(1-pval_px);    % threshold image at p-value

        [H,P,CI,STATS] = ttest2(signal_trials,noise_trials,'alpha',pval_px,'dim',1,'tail','right','vartype','unequal'); % one-sided, looking only for positive clusters
        real_tvals = STATS.tstat;
  
        %% METAPERMUTATIONS
        n_meta  = 20;
        n_permutes = 1000;
        
        ntrialsS = size(signal_trials,1);
        ntrialsN = size(noise_trials,1);
            
        cluster_thresh_bootc = zeros(1, n_meta);
        real_z_bootc = zeros(n_meta, 128*128);
        
        % centre data at each pixel for Wilcox's bootstrap-t - very important 
        csignal_trials = bsxfun(@minus,signal_trials,mean(signal_trials,1));
        cnoise_trials  = bsxfun(@minus,noise_trials,mean(noise_trials,1));
        
        tic
        for metapermi = 1:n_meta
            disp(metapermi);
            
            % initialize H0 matrices
            bootc_tvals    = zeros(n_permutes,128*128);      % Wilcox's bootstrap-t (resample with replacement within each condition): t-value distribution based on the real data)
            max_cluster_sizes_bootc = zeros(1,n_permutes);
            
            % generate pixel-specific H0 parameter distributions
            parfor permi = 1:n_permutes
                
                %% 3. BOOTSTRAPIN-t approach (mean-center each condition, sample with replacement, recompute t-statistic - allows to get data-based t-distribution)
                new_trialS_idx = randi(ntrialsS,ntrialsS,1);
                new_trialN_idx = randi(ntrialsN,ntrialsN,1);
                
                [H,P,CI,STATS] = ttest2(csignal_trials(new_trialS_idx,:),cnoise_trials(new_trialN_idx,:),'alpha',pval_px,'dim',1,'tail','right','vartype','unequal');
                bootc_tvals(permi,:) = STATS.tstat;
            end
           
            %**************************************************************
            % Wilcox's bootstrap-t
            %**************************************************************
            mean_h0_bootc = squeeze(nanmean(bootc_tvals));
            var_h0_bootc = squeeze(nanstd(bootc_tvals));
            
            % loop through permutations to find cluster sizes under the null hypothesis
            for permi = 1:n_permutes
                
                threshimg = squeeze(bootc_tvals(permi,:));
                threshimg = (threshimg-mean_h0_bootc)./var_h0_bootc;
                threshimg(isnan(threshimg))=0;
                
                % pixel-level statistics
                threshimg(threshimg<zval) = 0;        % only POSITIVE
                threshimg = reshape(threshimg, 128, 128,[]); % reshape
                
                % find clusters 
                islands = bwconncomp(threshimg);
                if numel(islands.PixelIdxList)>0
                    
                    % SUM
                    tempclustersum = zeros(1,length(islands.PixelIdxList));
                    for i = 1:numel(islands.PixelIdxList)
                        tempclustersum(i) = sum(abs(threshimg(islands.PixelIdxList{i})));
                    end
                    max_cluster_sizes_bootc(permi) = max(tempclustersum);
                end
            end
            %**************************************************************
            
            % find clusters and store it
            cluster_thresh_bootc(metapermi) = prctile(max_cluster_sizes_bootc,100-(100*pval_cluster)); % #RG: shoudl it be here divided too??
            
            % threshold real data
            real_z_bootc(metapermi,:) = (real_tvals-mean_h0_bootc)./var_h0_bootc;
        end
        toc
        
        %% get mean z-score from all metapermutations
        real_t_thresh_bootc = squeeze(mean(real_z_bootc,1));

        % threshold image at p-value
        real_t_thresh_bootc(real_t_thresh_bootc < zval) = 0;
        real_t_thresh_bootc(isnan(real_t_thresh_bootc)) = 0;             % turn NaNs to 0
        real_t_thresh_bootc = reshape(real_t_thresh_bootc, 128, 128,[]); % reshape

        %******************************************************************
        % Wilcox's bootstrap-t cluster-based testing
        %******************************************************************
        [islands,numclust] = bwlabel(real_t_thresh_bootc);
        
        for i=1:numclust
            % if cluster size is < than statistical threshold - remove it 
            tempcluster = sum(abs(real_t_thresh_bootc(islands(:)==i)));
            if tempcluster < mean(cluster_thresh_bootc)
                real_t_thresh_bootc(islands==i)=0;
            end
        end
        
        % ------------------------------------------------------------------
        % PLOT
        %------------------------------------------------------------------
        clf
        figure(1)
        
        alphadata = logical(real_t_thresh_bootc);
        clim = max(abs(mean(real_z_bootc,1)));
        
        im_main = imagesc(squeeze(reshape(real_tvals,128,128,[])));
        set(im_main,'AlphaData',alphadata);
        colormap(brewermap(256, '*RdYlBu'));
        hold on
        % Overlay Allen Map
        plot(Bregref(1),Bregref(2),'b*','markersize',5,'markerfacecolor','k')
        for p = 1:length(mymap.edgeOutline)-1
            plot(mymap.edgeOutline{p}(:, 2)-dx, mymap.edgeOutline{p}(:, 1)-dy,'k-');           % aligning Allen Bregma&Lambda with the Individual Bregma&Lambda
        end
        plot(mymap.bregma(1)-dx,mymap.bregma(2)-dy,'k*','markersize',5,'markerfacecolor','k')
        axis square
        axis square
        set(gca,'Ydir','reverse','xlim',[1 128],'ylim',[1 128],'clim',[-clim clim])
        set(gcf,'color','w')
        title([num2str(freqlist(freqi)) ' Hz, Wilcox, nperm = 1000, nmeta = 20'])
        colorbar
        
        figure_name = ([ writedir num2str(freqlist(freqi),'%02.f') '_Hz_boot_t' ]);
        saveas(gcf, figure_name, 'jpeg')
        %% Store the data
        
        real_t_thresh_bootc_all(freqi,:,:) = real_t_thresh_bootc;
        real_t_all(freqi,:,:) = reshape(real_tvals,128,128);
        
    end % END loop over frequencies
 
    %% Save 
    outfile = [writedir '\Tvals4clustering_boot_t.mat'];
    save(outfile,'real_t_all','real_t_thresh_bootc_all','pval_px','pval_cluster','zval','n_meta','n_permutes')

end % END loop over subjects


