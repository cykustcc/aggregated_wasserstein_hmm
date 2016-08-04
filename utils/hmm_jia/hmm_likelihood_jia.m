%  Copyright (c) 2015
%      Yukun Chen <cykustc@gmail.com>
% 
% Compute likelihood of given data sequence on given HMM model using Forward Backward Algorithm
%
% likelihood = hmm_likelihood_jia(data, d, nseq, onelen, numst, verbose, a00, a, gmm.mu, gmm.covariance, gmm.inv_covariance)
%
% Input:
%   - d: dimension
%   - nseq: number of sequence
%   - onelen: length for each sequence (Assume each sequence has the same length)
%	- numst: number of states
%   - verbose: if print out the estimated model 1: yes, 0: no
%   - a00: stationary distribution
%   - a: transition matrix
%   - gmm: gauss components
%       - gmm.mu: means of gauss components.
%       - gmm.covariance: covariances of guass components.
%       - gmm.inv_covariance: inv_covariances of guass components.
% Output:
%   - likelihood: the likelihood of the given data(sequence) on input gmmhmm model.
%
%
%
function likelihood = hmm_likelihood_jia(data, d, nseq, onelen, numst, verbose, a00, a, gmm.mu, gmm.covariance, gmm.inv_covariance)

    likelihood = hmm_likelihood_jia(data, d, nseq, onelen, numst, verbose, a00, a, gmm.mu, gmm.covariance, gmm.inv_covariance);
end