%% These scripts accompany the manuscript:
%  Gulbinaite et al. (2024) "Spatiotemporal resonance in mouse visual
%  cortex" Curr Biol
%
%  MATLAB version used 2019b
%% Code for spatiotemporal pattern analysis that employes two-stage GED 
% (described in Cohen (2022) NeuroImage PMID: 34906717 DOI:
% 10.1016/j.neuroimage.2021.118809)
 
% The two stages are: 1) dimensionality reduction using PCA;
%                     2) generalized eigenvalue decomposition (GED) or 
%                        Rhythmic Source Separation (RESS), see 
%                        Cohen and Gulbinaite, 2017 NeuroImage; PMID: 27916666)   

%  Requirements: 
%   - MATLAB Image Processing Toolbox; 
%   - Additional files in the folder: filterFGx.m

% INPUT: Preprocessed files (from 1_Preprocessing_singleTrial.m)  

% OUTPUT: Empirical spatial filters (spatiotemporal patterns) that  maximally 
% separate the signal data (response to specific flicker frequency) from 
% the reference data (brain response to spectrally adjacent non-flicker frequencies). 


%% PCA done only on the visual areas
clear all
close all

%% Directories 
rootdir = 'C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Data\'; 
subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

%% Load Allen CCF and mask
load('C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Code\mymap_30deg.mat')

allenmask = imread('C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Code\Allen_CCF_mask.tif','tiff');
allenmask = double(allenmask);

% Allen map info
Bregref = [10 30];
dx = mymap.bregma(1)-Bregref(1); % -5-Bregref(1);
dy = mymap.bregma(2)-Bregref(2);

%% Flicker frequencies used 

freqlist = [2:1:15 18:3:48 52:4:56 64:4:72];
%% ========================================================================
%                       Stage 1: PCA  
%% ========================================================================

%% LOOP over subjects
for subji = 1:length(subjects)
    disp(subji)

    cd([rootdir subjects(subji).name '\Ready4preprocessing\'])
    maskfile = dir('Mask*');
    mask0 = imread(maskfile.name,'tiff');
    mask0 = double(mask0);
     
    homedir = [rootdir subjects(subji).name '\Preprocessed_singleTrialGSR_RR_Gaussian\'];
    writedir = [rootdir subjects(subji).name '\PCA\'];
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end

    cd(homedir)

    %% LOOP over conditions (flicker frequencies)
    for freqi = 1:length(freqlist)
        clear subdata tmpdata V D

        %% Select files with appropriate flicker freq
        
        filelist = dir(['*' num2str(freqlist(freqi),'%02.f') '_Hz*']);
        filelist = {filelist.name};
        
        %% Outputfile
        
        outfilename = [ writedir 'PCA' filelist{1}(11:16) '.mat'];
        if exist(outfilename,'file'), continue; end
        
        %% Load data
       
        for filei = 1:length(filelist)
            disp(filelist{filei})
            filedata = load([ homedir filelist{filei}]);
  
            if filei == 1
                data2analyze = zeros(length(filelist),128*128,length(filedata.time4erp));
            end
            
            % Reapply the mask, because Gaussian filter "extends" the masked area 
            filedata.trialdata_3D = bsxfun(@times, filedata.trialdata_3D, mask0);
          
            % Store in 2D
            data2analyze(filei,:,:) = reshape(filedata.trialdata_3D,128*128,[]);
        end
    
       NonZeroPixelsIndex = filedata.NonZeroPixelsIndex;
       %% Take only flicker-part of the data
       
       num_px = 128*128;       
       time4ssvep = filedata.time4erp;
       time4fft  = [1000 10000]; % skip the first 1000 ms due to transient response; same as in 2_FFT_power_angle_maps.m
       tidx      = dsearchn(time4ssvep',time4fft');
       
       data2use = data2analyze(:,:,tidx(1):tidx(2));
       
       %% Change dimensions to pixels x time x trials
       
       num_px = 128*128;
       num_pnts = size(data2use,3);
       ntrials = size(data2use,1);
       
       tmpdata = zeros(num_px,num_pnts,ntrials);
       parfor pixi = 1:num_px
           tmpdata(pixi,:,:) = squeeze(data2use(:,pixi,:))';
       end
       
       %% PCA for generic dimensionality reduction.
       
       tmpdata = tmpdata(NonZeroPixelsIndex,:,:);
      
       % covariance from broadband data
       covA = zeros(length(NonZeroPixelsIndex));
       for triali=1:ntrials
           tmp = detrend(  tmpdata(:,:,triali)' )';
           covA = covA + tmp*tmp'/num_pnts;
       end
       
       % Perform PCA and sort eigen vectors 
       tic
       [V,D]  = eig(covA/ntrials);
       toc
       [D,sx] = sort(diag(D),'descend');
       V = V(:,sx);     % sort eigvecs
       
       % convert evals to %change
       D = 100*D/sum(D);
       
       % which components to keep (threshold of % total variance accounted for)
       comps2keep = D>.01; %was:0.01
       nComps = sum(comps2keep);
       
       % reconstruct compressed data: components x time x trials
       subdata = reshape( (reshape(tmpdata,length(NonZeroPixelsIndex),[])' * V(:,comps2keep))' ,[nComps num_pnts ntrials]);
       subdata = single(subdata);
       V = single(V);
       
       %% save PCA component data
       save(outfilename,'subdata','nComps','NonZeroPixelsIndex','num_pnts','ntrials','num_px','V','time4fft','time4ssvep')

    end
end

%% ========================================================================
%             Stage 2: GED or Rhythmic Source Separation (RESS)
%% ========================================================================

set(0,'DefaultFigureWindowStyle','normal')

clear all
close all

% Directories
rootdir = 'C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Data\'; 
subjects = dir(rootdir);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

%% Load Allen CCF and mask
load('C:\Users\abete\Desktop\Spatiotemporal_resonance_DRYAD\Code\mymap_30deg.mat')

% Allen map info
Bregref = [10 30];

dx = mymap.bregma(1)-Bregref(1); % -5-Bregref(1);
dy = mymap.bregma(2)-Bregref(2);

%% Flicker frequencies used 

freqlist = [2:1:15 18:3:48 52:4:56 64:4:72];

%% LOOP over subjects
for subji = 1%:length(subjects)
    disp(subji)

    %% Directories
    
    cd([rootdir subjects(subji).name '\Ready4preprocessing\'])
    maskfile = dir('Mask*');
    mask0 = imread(maskfile.name,'tiff');
    mask0 = double(mask0);
    
    homedir = [rootdir subjects(subji).name '\PCA\'];
    writedir = [rootdir subjects(subji).name '\RESS\'];
    if(exist(writedir, 'dir') == 0), mkdir(writedir); end

    cd(homedir)

    %% LOOP over frequencies
    for freqi = 5%1:length(freqlist)
        
        %% Load data after PCA
        
        filename = dir(['*' num2str(freqlist(freqi),'%02.f') '_Hz*']);
        disp(filename.name)
        load([ homedir filename.name]);

        %% Initialize matrices to store the results
        
        if freqi== 1
            resstopos = zeros(length(freqlist),10,128*128);            
            ress = zeros(length(freqlist),num_pnts,ntrials);
            resonance_profile = zeros(1,length(freqlist));
            signal_noise_SNR  = zeros(length(freqlist),length(freqlist));
            
            resonance_evals = zeros(1,length(freqlist));
            var_explained = cell(1,length(freqlist));  % evals expressed as %change
        end
        
        %% Temporal filter around flicker and neighbouring frequencies
        
        srate = 150.64;
        freq2use = freqlist(freqi);
        peakwidt  = 0.5; % in Hz
        neighfreq = 2;   % in Hz
        neighwidt = 2;   % in Hz
        tidx = [1 num_pnts];
        
        % get covariance for flicker frequency
        filtTmp = filterFGx(subdata,srate,freq2use,peakwidt);
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),nComps,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covAt = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        % get covariance for below flicker frequency
        filtTmp = filterFGx(subdata,srate,freq2use-neighfreq,neighwidt);
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),nComps,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covLo = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        % get covariance for above flicker frequency
        filtTmp = filterFGx(subdata,srate,freq2use+neighfreq,neighwidt);
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),nComps,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covHi = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% GED
        
        covR = (covHi+covLo)/2;
        
        %% Shrinkage regularization of R covariance matrix

        gamma = .01;
        covRr = covR*(1-gamma) + eye(nComps)*gamma*mean(eig(covR));
        [evecs,evals] = eig(covAt,covRr);

        %% sort eigenvectors/eigenvalues
        [evals, sidx] = sort(diag(evals), 'descend');
        evecs = real(evecs(:, sidx));
        
        evecs_unnormed = evecs;
        evecs_rel = zeros(size(evecs));
        for v = 1:size(evecs,2)
            evecs(:,v) = evecs(:,v)/norm(evecs(:,v)); % normalize to unit length
            rel_eval = evals(v)/sum(evals);           % extract relative eigenvalue
            evecs_rel(:,v) = evecs(:,v) * rel_eval;   % normalize to relative eigenvalue
        end
        
        % convert evals to %change
        evals_percent = 100*evals/sum(evals);

        % PLOT
        %---------------
        f = figure('Position', get(0, 'Screensize'));
        set(gcf,'name', 'neigh = 2Hz; FWHM = 2Hz')
        
        subplot(271)
        plot(evals,'-o')
        title( ['Eigenvalues: ' num2str(freq2use) 'Hz' ])
        axis square

        %% Reconstruct RESS timeseries (10 first components)

        nRESS = 10;
        tmpresstx = zeros(nRESS,num_pnts,ntrials);
        topos = zeros(nRESS,num_px);
        for triali=1:ntrials
            tmpresstx(:,:,triali) = ( subdata(:,:,triali)' * evecs(:,1:nRESS) )'; 
        end
        
        % Compute RESS topomaps: note the additional transformation using PCA eigenvectors
        V = double(V);
        for compi = 1:nRESS
             topos(compi,NonZeroPixelsIndex) = real(evecs(:,compi)'*covAt*V(:,1:nComps)');  
        end

        % fix the sign of components
        for ci=1:nRESS
            [~,idx] = max(abs(topos(ci,:)));                 % find component with highest weight
            topos(ci,:) = topos(ci,:) * sign(topos(ci,idx)); % force to positive sign
        end

        % Compute RESS component power spectra
        resolution = 0.1;
        nFFT = ceil( srate/resolution );
        hz   = srate/2*linspace(0,1,floor(nFFT/2+1));

        skipbins = round(0.5/diff(hz(1:2))); % skip 0.5 Hz
        numbins  = round(1/diff(hz(1:2))) + skipbins; % +-1 Hz = 10
        
        freqidx = dsearchn(hz',freq2use);
        
        tmpresspwr = mean(abs(fft(tmpresstx,nFFT,2)).^2,3);
        tmpresspwr = squeeze(tmpresspwr(:,floor(1:nFFT/2+1)));
        
        % Translate pow to SNR
        for compi=1:nRESS
            for hzi=numbins+1:length(hz)-numbins-1
                numer = squeeze(tmpresspwr(compi,hzi));
                denom = mean( tmpresspwr(compi,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]),2 );
                tmpresspwr(compi,hzi) = numer./denom;
            end
        end
        
        % Store power at all flicker frequencies used - used later for noise SNR
        % (H0 or SNR expected when flicker of harmonic frequenies were not presented
        for i = 1:length(freqlist)
            freq2use_extra = freqlist(i);
            freqidx_extra = dsearchn(hz',freq2use_extra);
            
            signal_noise_SNR(freqi,i)  = squeeze(tmpresspwr(1,freqidx_extra)); % take only the 1st component
        end
        
        % Store SNRs from each component
        tmpsnr = squeeze(tmpresspwr(:,freqidx));

        % Plot the first 6 RESS components
        for compi = 1:6
            
            subplot(2,7,compi+1)
            imagesc(reshape(squeeze(topos(compi,:)),128,128,[]));
            hold on
            % Overlay Allen Map
            plot(Bregref(1),Bregref(2),'b*','markersize',5,'markerfacecolor','w')
            for p = 1:length(mymap.edgeOutline)-1
                plot(mymap.edgeOutline{p}(:, 2)-dx, mymap.edgeOutline{p}(:, 1)-dy,'w-');           
            end
            plot(mymap.bregma(1)-dx,mymap.bregma(2)-dy,'w*','markersize',5,'markerfacecolor','w')
            set(gca,'xlim',[20 75],'ylim',[55 110])
            axis square
            title( ['RESS: ' num2str(compi)])

            ax(compi+8) = subplot(2,7,compi+8);
            plot(hz,squeeze(tmpresspwr(compi,:)) ,'b','Linewidth',1.5);
            title( ['SNR = ' num2str(round(tmpsnr(compi))) ])

            axis square
            set(gca,'xlim',[freq2use-4 freq2use+4],'ylim',[0 tmpsnr(compi)*1.1])
        end
        set(gcf,'color','w')
        
        % save figure
        figure_name = (['RESS_topos_' filename.name(5:9)]);
        saveas(gcf, [writedir  figure_name], 'jpeg')

        %% save RESS results
        comp2plot = 1;
        resonance_profile(freqi) = tmpsnr(comp2plot);
        resonance_evals(freqi)   = evals(comp2plot);
        var_explained{freqi}     = evals_percent;
        resstopos(freqi,:,:)     = topos;
        
        if subji == 6 && freqi == 18 % the last trial of subj 6 Rec 5 is missing post-stimulus period and thus was not included
            ress(freqi,:,1:9) = squeeze(tmpresstx(1,:,:));
        else
            ress(freqi,:,:) = squeeze(tmpresstx(1,:,:));
        end

    end % END loop over frequencies
 
    %% SAVE results
    outfilename = [writedir subjects(subji).name '_RESS.mat'];
    resstopos = single(resstopos);
    ress = single(ress);
    save(outfilename,'var_explained','resstopos','signal_noise_SNR','ress','resonance_profile','resonance_evals','hz','nFFT','NonZeroPixelsIndex','peakwidt','neighfreq','neighwidt','resolution','skipbins','numbins')
    
    %% PLOT and save figures
    figure
    subplot(121)
    plot(freqlist,resonance_profile,'ro-','Linewidth',2,'MarkerFace','w')
    set(gca,'xlim',[1.8 79],'xscale','log', 'xtick',round(logspace(log10(freqlist(1)),log10(freqlist(end)),6)))
    legend({ 'Comp = 1'})
    ylabel('Power (SNR)')
    xlabel('Flicker frequencies (Hz)')
    axis square
    box off
    legend box off
    title('Renonace profile (RESS)')
    
    subplot(122)
    plot(freqlist,resonance_evals,'ro-','Linewidth',2,'MarkerFace','w')
    set(gca,'xlim',[1.8 79],'xscale','log', 'xtick',round(logspace(log10(freqlist(1)),log10(freqlist(end)),6)))
    set(gcf,'color','w')
    axis square
    title('Max eigenvalue profile')
    xlabel('Flicker frequencies (Hz)')
    ylabel('Eigenvalue')
    box off
    set(gcf,'color','w')
    
    figure_name = ['Resonance_profile_RESS_' subjects(subji).name];
    saveas(gcf, [writedir figure_name], 'fig')

end % END loop over subjects
