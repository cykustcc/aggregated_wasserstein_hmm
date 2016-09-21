function dist = wasserstein_gaussian( mu1,cov1,mu2,cov2 )
%WASSERSTEIN_GAUSSIAN Returns the norm 2 wasserstein distance of 2
%gaussians
% Authors: Yukun Chen, College of Information Science and Technology, Penn State University
% Contact: cykustc@gmail.com
%    cov1=(cov1+cov1')/2;
%    cov2=(cov2+cov2')/2;
%     w2=norm(mu1-mu2)^2+trace(cov1+cov2-2*sqrtm(sqrtm(cov1)*cov2*sqrtm(cov1)));
    sqrt_cov1=sqrtm(cov1);
    w2=norm(mu1-mu2)^2+trace(cov1+cov2-2*sqrtm(sqrt_cov1*cov2*sqrt_cov1));
    dist=real(w2^0.5);
%     dist=real(w2);
end

