function R = GSR_NonZeroPixels(P)
% function y = GSR(stack)
%% Global signal regression

% Edgar Bermudez October 2016
% adapted from Fox MD, Zhang D, Snyder AZ, Raichle ME. The global signal
% and observed anticorrelated resting state brain networks. J Neurophysiol
% [Internet]. 2009;101(6):3270–83. 

% Edited by: Rasa Gulbinaite March 2021, adapted for 2D instead of 3D matrix

% stack is a pixels x frames matrix
% height = size(stack,1);     % #RG: commented out
% width = size(stack,2);      % #RG: commented out
% num_frames = size(stack,3); % #RG: commented out

num_frames = size(P,2);       % #RG: P = pixels x frames

% P = reshape(stack, height*width, num_frames); % #RG no need to reshape, P is a 2D matrix 

% global signal regression
% Pr = g^+ * P
% the regressed pixels Pr are the product of the pseudo-inverse g^+ times
% the pixel x frames matrix P. The global mean time course, g, is a
% pixel x 1 vector, g = (1/num_frames)B * 1m where 1m is the m x 1 vector
% of all ones.

% the regresion g^+ = (g' * g)^-1 * g' and g^+ * g = 1;

% calculate g
g = (1/num_frames)*(P*ones(num_frames, 1));

% g^+
g_p = (g'*g)\g';

% regression for Pixels P using g
P_g =  g_p * P;

% regressed pixels 
R = P - (g*P_g);


% y = reshape(R, height, width, num_frames); % #RG: do not need to reshape back in 2D version 

