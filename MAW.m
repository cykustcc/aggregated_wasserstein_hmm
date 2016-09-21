function dist = MAW( gmmhmm1, gmmhmm2, alpha)
% Compute Minimized Aggregated Wasserstein (MAW) distance between hidden 
% Markov Models with Gaussian Mixture emission functions (GMMHMM) with 
% different methods
% Input:
% gmmhmm1, gmmhmm2 -- objects that stores the parameters for the two 
%                     GMMHMMs in comparison. (See gmmhmm.m for how to 
%                     construct such an object with GMMHMM parameters)
% alpha            -- parameter for weighted sum of the difference between
%                     two marginal GMMs of two GMMHMMs and the difference 
%                     between the transition matrices of two GMMHMMs.
% option           -- shared:
%                       option.method: {'badmm', 'sinkhorn'}
%
%                     method specific options:
%
    D = pdist2(gmmhmm1, gmmhmm2, 'gaussianMixture');
    [dist, matching]=gmm_MAW(gmmhmm1,gmmhmm2, D);
    dist2 = gmmhmm_MAW_transmat(gmmhmm1,gmmhmm2, matching, D);
    dist = (1-alpha)*dist + alpha*dist2;
end