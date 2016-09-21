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
%                     IAW = (1-alpha)*gmm_diff + alpha*trans_mat_diff
% optoins          -- shared:
%                       options.method: {'badmm', 'sinkhorn'}
%                       options.sample_size: size of samples generated from
%                                           two marginal GMMs, (default:
%                                           100)
%                       options.max_iter: max iteration for solving the 
%                                        optimal transport problem. (defa-
%                                        -ult: 300 )
%
%                     method specific options:
%                       'sinkhorn':
%                       options.sinkhorn_epsilon (default: 0.01)
%                       options.sinkhorn_max_iters (default: 300)
%                       'badmm':
%                       options.badmm_max_iters (default: 300)


    method = 'sinkhorn';
    if isfield(options, 'method')
        method = options.method;
    end
    
    sample_size = 100;
    if isfield(options, 'sample_size')
       sample_size = options.sample_size; 
    end
    
    if strcmp(method, 'sinkhorn') 
        [dist, matching] = gmm_IAW_Sinkhorn(gmmhmm1,gmmhmm2,sample_size, options);
    elseif strcmp(method, 'badmm')
        [dist, matching] = gmm_IAW_BADMM(gmmhmm1,gmmhmm2,sample_size, options);
    end 
    
    if any(isnan(matching(:)))
        disp('matching has nan.')
        error('%f',matching);
    end
    D = pdist2(gmmhmm1, gmmhmm2,'gaussianMixture');
    [dist2]=gmmhmm_MAW_transmat(gmmhmm1,gmmhmm2,matching, D);
    dist = dist + dist2;
end