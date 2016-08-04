%  Copyright (c) 2015
%      Yukun Chen <cykustc@gmail.com>
% 
% Fit HMM model (Gaussian emission function) using E-M algorithm
%
% [a00, a, gauss] = hmm_fit_jia(V, N, source)
%
% Input:
%   - d: dimension
%   - nseq: number of sequence
%   - onelen: length for each sequence (Assume each sequence has the same length)
%	- numst: number of states
%   - verbose: if print out the estimated model 1: yes, 0: no
%
% Output:
%   - a00: stationary distribution
%   - a: transition matrix
%   - guass: gauss components
%       - gauss.means: means of gauss components.
%       - gauss.sigma: covariances of guass components.
%
% The distance is approximated by the graph distance between points.
% The path and distance are computed with Dijkstra's algorithm.
%
%
function [a00, a, gauss] = hmm_fit_jia(data, d, nseq, onelen, numst, verbose)

    [a00, a, gauss] = hmm_fit_jia(data, d, nseq, onelen, numst, verbose);
end