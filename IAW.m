function dist = IAW( gmmhmm1, gmmhmm2, alpha, options)
% Compute Improved Aggregated Wasserstein (IAW) distance between hidden 
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
%                       option.sample_size: size of samples generated from
%                                           two marginal GMMs, (default:
%                                           100)
%
%                     method specific options:
%
    method = 'sinkhorn';
    if isfield(options, 'method')
        method = options.method;
    end
    
    sample_size = 100;
    if isfield(options, 'sample_size')
       sample_size = options.sample_size; 
    end
    
    if strcmp(method, 'sinkhorn') 
        [dist, matching] = gmm_IAW_Sinkhorn(gmmhmm1,gmmhmm2,sample_size);
    elseif strcmp(method, 'badmm')
        [dist, matching] = gmm_IAW_badmm(gmmhmm1,gmmhmm2,sample_size);
    end 
    
    if any(isnan(matching(:)))
        disp('matching has nan.')
        error('%f',matching);
    end
    D = pdist2(gmmhmm1, gmmhmm2,'gaussianMixture');
    [dist2]=gmmhmm_MAW_transmat(gmmhmm1,gmmhmm2,matching, D);
    dist = dist + dist2;
end