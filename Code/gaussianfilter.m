function y = gaussianfilter(data, hsize, sigma)
%% gaussianfilter(data, hsize, sigma) to smooth in 2D 
% Edgar Bermudez
% edgar.bermudez@uleth.ca


gaus_filt = fspecial('gauss', hsize, sigma);

height = size(data,1);
width = size(data,2);
num_frames = size(data, 3);

filtdata = zeros(size(data));

parfor fr = 1:num_frames
	filtdata(:,:,fr) = conv2(data(:,:,fr), gaus_filt, 'same');
end

y = filtdata;
