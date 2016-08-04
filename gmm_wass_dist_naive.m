function [fval, x, lambda] = gmm_wass_dist_naive( gmm1,gmm2, D )
%GMM_WASS_NAIVE Summary of this function goes here
%   Detailed explanation goes here
% [f,PI,lambda]=wasserstein(gmm1.weights,gmm2.weights);

%     optim_options   = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off');
%     lpoptim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off', 'Simplex', 'on');
%     default_options = optimset('Display','off', 'Diagnostics','off');
%    global optim_options lpoptim_options default_options;

    wX=gmm1.weights;
    wY=gmm2.weights;

    wX=wX+eps;
    wY=wY+eps;

    wX = wX/sum(wX);
    wY = wY/sum(wY);
    %if any(wX<eps) || any(wY<eps)
    %    error('%f ',[wX wY]);
    %end

%     D = pdist2(X', Y', 'seuclidean');
    if nargin <3
    D = pdist2(gmm1,gmm2,'gaussianMixture');
    end

%     [x, fval, exitflag, ~, lambda] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], x0, default_options );
%    [x, fval, exitflag, ~, lambda] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], x0, lpoptim_options );

    [exitflag, fval, x, lambda] = OptimalTransport_LP(D, wX', wY');

    x = reshape(x, size(D));
    x(x<0) = 0;
end

