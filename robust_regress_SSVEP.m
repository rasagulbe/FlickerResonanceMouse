function robust_regress_SSVEP(triali, filename, outfilename)

fprintf([ '\n\nStarting trial ' num2str(triali) '.' ]);

%% Check if outputfile exists

if exist(outfilename,'file'), return; end
%% load data

load(filename,'trialdata','NonZeroPixelsIndex','time4erp'); 
%% Robust regression settings
basis = 'polynomials';
thresh = 3;
niter = 3;
wsize = 100;
order = 20;

camera_srate = 150;
start_mask   = round(0.0*camera_srate); % in samples/frames 
end_mask     = round(0.7*camera_srate); % in samples/frames
start_mask2  = round(10*camera_srate);
end_mask2    = round(10.7*camera_srate);

% Create a temporal mask for the trial (same for all trials)
trial_temporal_mask2 = ones(1,size(trialdata,2));  % array to store the temporal mask to mask-out all trial onsets

zerodix = dsearchn(time4erp',0);
mask_startidx  = zerodix - start_mask;
mask_stopidx   = zerodix + end_mask;
mask2_startidx = zerodix + start_mask2; 
mask2_stopidx  = zerodix + end_mask2;

trial_temporal_mask2(mask_startidx:mask_stopidx) = 0;
trial_temporal_mask2(mask2_startidx:mask2_stopidx) = 0;
%% Robust regression on NonZeropixels

% weights constrain which parts of the signal will be fitted; in non-fitted or masked parts w = 0
wt = repmat(trial_temporal_mask2,[length(NonZeroPixelsIndex) 1]);

% call nt_detrend with ERP part masked out by wt
[tmp,w1,regressline1] = nt_detrend_edited(trialdata',1,wt',basis,thresh,niter,wsize); % start with 1st order

% subtract the regression line
trialdata = trialdata - regressline1';

% call nt_detrend with higher order polynomial
[yy,ww,regressline2] = nt_detrend_edited(tmp,order,w1,basis,thresh,niter,wsize); % then nth order with mask of previous step

% subtract the regression line
trialdata = trialdata-regressline2';

%% Save data
save(outfilename,'trialdata','wsize','order','NonZeroPixelsIndex','time4erp','start_mask','start_mask2','end_mask','end_mask2')

end