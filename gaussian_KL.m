function dist = gaussian_KL( mu1,cov1,mu2,cov2 )
%gaussians
% Authors: Yukun Chen, College of Information Science and Technology, Penn State University
% Contact: cykustc@gmail.com
%     dist=1.0/2.0*((mu1-mu2)'*cov1*(mu1-mu2)+trace(inv(cov1)*cov2)-length(mu1)+log(det(cov1)/det(cov2)));
    cov1=(cov1+cov1')/2;
    cov2=(cov2+cov2')/2;
    inv_cov2=inv(cov2);
    dist=1.0/2.0*((mu1-mu2)'*inv_cov2*(mu1-mu2)+trace(inv_cov2*cov1)-length(mu1)+log(det(cov2)/det(cov1)));
%     dist=log(det(cov2)/det(cov1))+(det(cov1)^2+norm((mu1-mu2))^2)/det(cov2)^2/2;
end