function [ distance ] = gmmhmm_KL( gmmhmm1,gmmhmm2,sample_size )
% gmmhmm likelihood distance, 
% The query sample should be gmmhmm2
% for more details, please refer to: B.-H. F. Juang and L. R. Rabiner, A probabilistic distance measure for hidden markov models,? AT&T technical journal, vol. 64, no. 2, pp. 391?408, 1985. 
data=gmmhmm1.gen_mk_chain(sample_size);
numst = gmmhmm1.mode_num;
distance = hmm_likelihood_jia(data, gmmhmm1.dim, 1, sample_size,numst,0,gmmhmm1.weights, gmmhmm1.trans_mat, gmmhmm1.mu, gmmhmm1.covariance, gmmhmm1.inv_covariance);
distance = distance - hmm_likelihood_jia(data, gmmhmm2.dim, 1, sample_size,numst,0,gmmhmm2.weights, gmmhmm2.trans_mat, gmmhmm2.mu, gmmhmm2.covariance, gmmhmm2.inv_covariance);
distance = distance / size(data,2);
end

