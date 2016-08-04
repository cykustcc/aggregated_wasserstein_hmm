mkdir_if_not_exist([gmmhmm_projectroot,'/imgs/wass_vs_KL/']);
SAMPLE_SIZE=50;
SEQ_NUM=100;
cov_factor=1;
eg.dim=2;
eg.mu=[[0,0]]';
eg.covariance=cat(3,cov_factor*[1,0;0,1]);
gmm_eg0 = gmm(eg.dim,1,eg.mu,eg.covariance,[1]);
a=gmm_eg0.rndlist(SAMPLE_SIZE);

gmm_eg0.plot_rndsamples(50);
print([gmmhmm_projectroot,'/imgs/wass_vs_KL/','samples.png'], '-dpng','-r100');

% disp(a);
disp_factor=1;

GMM_NUM=10;
gmm_model=[];
wass_Dist=zeros(SEQ_NUM,GMM_NUM);
KL_Dist=zeros(SEQ_NUM,GMM_NUM);
for i=1:GMM_NUM
    egi.dim=2;
    egi.mu=[[i*disp_factor,i*disp_factor]]';
    egi.covariance=cat(3,cov_factor*[1,0;0,1]);
    gmm_egi = gmm(egi.dim,1,egi.mu,egi.covariance,[1]);
    for j=1:SEQ_NUM
        data=gmm_eg0.rndlist(SAMPLE_SIZE);
        obj=gmdistribution.fit(data',1);
        estimated_gmm_eg0=gmm(eg.dim,1,obj.mu',cat(3,obj.Sigma),[1]);
 
        wass_Dist(j,i)=gmm_wass_dist_naive(gmm_egi,estimated_gmm_eg0);
%         gt_wass_dist=gmm_wass_dist_naive(gmm_egi,gmm_eg0);
        
        
        KL_Dist(j,i)=gaussian_KL(gmm_egi.mu,gmm_egi.covariance,estimated_gmm_eg0.mu,estimated_gmm_eg0.covariance);
%         gt_KL_dist=gaussian_KL(gmm_egi.mu,gmm_egi.covariance,gmm_eg0.mu,gmm_eg0.covariance);
%         KL_Dist(j,i)=KL_Dist(j,i)^0.5*gt_wass_dist/gt_KL_dist^0.5*det(gmm_eg0.covariance);
        KL_Dist(j,i)=KL_Dist(j,i)^0.5*det(gmm_eg0.covariance);
    end
end
mkdir_if_not_exist([gmmhmm_projectroot,'/data/wass_vs_KL/']);
save([gmmhmm_projectroot,'/data/wass_vs_KL/','wass_Dist.mat'],'wass_Dist');
save([gmmhmm_projectroot,'/data/wass_vs_KL/','KL_Dist.mat'],'KL_Dist');
load([gmmhmm_projectroot,'/data/wass_vs_KL/','wass_Dist.mat']);
load([gmmhmm_projectroot,'/data/wass_vs_KL/','KL_Dist.mat']);

%% Compute distance and its error bar
wass_Dist_mean=mean(wass_Dist,1);
wass_Dist_var=3*std(wass_Dist,1);

KL_Dist_mean=mean(KL_Dist,1);
KL_Dist_var=3*std(KL_Dist,1);

%% Plot figure
figure;
hold on;
e1=errorbar(wass_Dist_mean,wass_Dist_var);
set(e1,'Linewidth',3)
hold on;
e2=errorbar(KL_Dist_mean,KL_Dist_var);
set(e2,'Linewidth',3)
grid on;
title(gca,['Distance value w.r.t deviation factor of \mu'])
xlabel('deviation factor i', 'fontsize', 25);
ylabel('Distance Value', 'fontsize', 25);
set(gca, 'linewidth', 3, 'fontsize', 20);
legend({'Wasserstein Distance','KL Divergence'}, 'location', 'northwest');
print([gmmhmm_projectroot,'/imgs/wass_vs_KL/','wass_vs_KL.png'], '-dpng','-r100');
print([gmmhmm_projectroot,'/imgs/wass_vs_KL/','wass_vs_KL.eps'], '-depsc','-r100');

