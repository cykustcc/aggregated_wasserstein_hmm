function [ distance ] = gmmhmm_dist_likelihood_raw( gmmhmm1,data )
% gmmhmm likelihood distance using raw data, 
% The query sample should be gmmhmm1
% for more details, please refer to: B.-H. F. Juang and L. R. Rabiner, A probabilistic distance measure for hidden markov models,? AT&T technical journal, vol. 64, no. 2, pp. 391?408, 1985. 
numst = gmmhmm1.mode_num;
% distance = hmm_likelihood_jia(data, gmmhmm_data.dim, 1, size(data,2),numst,0,gmmhmm_data.weights, gmmhmm_data.trans_mat, gmmhmm_data.mu, gmmhmm_data.covariance, gmmhmm_data.inv_covariance);
distance = -hmm_likelihood_jia(data, gmmhmm1.dim, 1, size(data,2),numst,0,gmmhmm1.weights, gmmhmm1.trans_mat, gmmhmm1.mu, gmmhmm1.covariance, gmmhmm1.inv_covariance);
distance = distance / size(data,2);
end
