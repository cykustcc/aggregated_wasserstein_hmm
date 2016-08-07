function fval = gmm_MAW_givenMatching( gmm1,gmm2,matching)
%GMM_WASS_NAIVE Summary of this function goes here
%   Detailed explanation goes here
% [f,PI,lambda]=wasserstein(gmm1.weights,gmm2.weights);

%     optim_options   = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off');
%     lpoptim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off', 'Simplex', 'on');
%     default_options = optimset('Display','off', 'Diagnostics','off');


%     D = pdist2(X', Y', 'seuclidean');
    D = pdist2(gmm1,gmm2,'gaussianMixture');
    
    dist=D.*matching;
    fval=sum(dist(:));
end



