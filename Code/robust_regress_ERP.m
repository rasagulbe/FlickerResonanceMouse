function robust_regress_ERP(triali, filename, outfilename,trial_temporal_mask,time2analyze)

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
order = 50; 

%% Other settings
time2analyze = [-2000 4000]; % smaller window for efficiency
time2analyze_idx = dsearchn(time4erp',time2analyze');
trialdata = trialdata(:,time2analyze_idx(1):time2analyze_idx(2));
time4erp = time4erp(time2analyze_idx(1):time2analyze_idx(2));

%% Robust regression on NonZeropixels

% weights define which parts of the signal will be fitted; masked parts have w=0;
wt = repmat(trial_temporal_mask,[length(NonZeroPixelsIndex) 1]);

% call nt_detrend with brain data masked out by wt
[tmp,w1,regressline1] = nt_detrend_edited(trialdata',1,wt',basis,thresh,niter,wsize); % start with 1st order

% subtract the regression line
trialdata = trialdata - regressline1';

% call nt_detrend with higher order polynomial
[yy,ww,regressline2] = nt_detrend_edited(tmp,order,w1,basis,thresh,niter,wsize); % then nth order with mask of previous step

% subtract the regression line
trialdata = trialdata-regressline2';

%% Save data
save(outfilename,'trialdata','wsize','order','NonZeroPixelsIndex','time4erp')

end
