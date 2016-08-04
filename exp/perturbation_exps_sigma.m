%%%%%%%%%%%%%%%%%%%
%% Perturb sigma
%%%%%%%%%%%%%%%%%%%
%% gmm_hmm_naive
factors = [0.2,0.4,0.6];
paras={};
tic
for factor_idx = 1:length(factors)
    para={};
    factor=factors(factor_idx);
    disp(factor);
    GMM_NUM=5;
    SEQ_NUM_FOR_EACH_GMM=10;
    gmmhmm_models=[];
    for i=1:GMM_NUM
        eg.dim=2;
        eg.mu=[[2,2];[5,5]]';
        eg.trans_mat=[.8 .2; .2 .8;];
        eg.weights=station_dist(eg.trans_mat);
        r_num=rand(2);
        r_num=(r_num+r_num')/2;
        temp = 0.2*expm(factor*i*r_num);
        eg.covariance=cat(3,temp,temp);
        gmmhmm_eg1 = gmmhmm(eg.dim,eg.trans_mat,eg.weights,eg.mu,eg.covariance);
        n=100; % length of chain
        for j=1:SEQ_NUM_FOR_EACH_GMM
            rng(i*j); % make result reproducible
            data = gmmhmm_eg1.gen_mk_chain(n);
            numst = size(eg.mu,2);
            [a_00, a, gauss]=hmm_fit_jia(data, eg.dim, 1, n,numst,0);
            means=vertcat(gauss.means)';
            cov=cat(3,gauss.sigma);
            est_a = gmmhmm(eg.dim,a,a_00,means,cov);
            gmmhmm_models=[gmmhmm_models,est_a];
        end
    end
    dist_matrix_dim=GMM_NUM*SEQ_NUM_FOR_EACH_GMM;
    dist_matrix=zeros(dist_matrix_dim);
    dist_matrix2=zeros(dist_matrix_dim);
    for i=1:dist_matrix_dim
        for j=1:i
            D = pdist2(gmmhmm_models(i),gmmhmm_models(j),'gaussianMixture');
            [dist,matching]=gmm_wass_dist_naive(gmmhmm_models(i),gmmhmm_models(j),D);
            dist2 = gmmhmm_dist_naive_transmatdist(gmmhmm_models(i),gmmhmm_models(j),matching,D);
            dist_matrix(i,j)=dist;
            dist_matrix(j,i)=dist;
            dist_matrix2(i,j)=dist2;
            dist_matrix2(j,i)=dist2;
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix,dist_matrix2,ground_truth_class,['MAW \Delta\sigma= ',num2str(factor)]};
    paras{factor_idx}=para;
end
toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
% paras{1}{4}='MAW \Delta\sigma= 0.2';
% paras{2}{4}='MAW \Delta\sigma= 0.4';
% paras{3}{4}='MAW \Delta\sigma= 0.6';
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_naive.mat'],'paras');
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_naive.mat'])
% draw_multiple_precrecl(paras,'Changing \Sigma','change_sigma_naive','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/',0.1)


%% KL divergence likelihood --Juang's paper
factors = [0.2,0.4,0.6];
paras={};
tic
for factor_idx = 1:length(factors)
    para={};
    factor=factors(factor_idx);
    disp(factor);
    GMM_NUM=5;
    SEQ_NUM_FOR_EACH_GMM=10;
    gmmhmm_models=[];
    for i=1:GMM_NUM
        eg.dim=2;
        eg.mu=[[2+i*factor,2+i*factor];[5+i*factor,5+i*factor]]';
        eg.trans_mat=[.8 .2; .2 .8;];
        eg.weights=station_dist(eg.trans_mat);
        r_num=rand(2);
        r_num=(r_num+r_num')/2;
        temp = 0.2*expm(factor*i*r_num);
        eg.covariance=cat(3,temp,temp);
        gmmhmm_eg1 = gmmhmm(eg.dim,eg.trans_mat,eg.weights,eg.mu,eg.covariance);
        n=100; % length of chain
        for j=1:SEQ_NUM_FOR_EACH_GMM
            rng(i*j); % make result reproducible
            data = gmmhmm_eg1.gen_mk_chain(n);
            numst = size(eg.mu,2);
            [a_00, a, gauss]=hmm_fit_jia(data, eg.dim, 1, n,numst,0);
            means=vertcat(gauss.means)';
            cov=cat(3,gauss.sigma);
            est_a = gmmhmm(eg.dim,a,a_00,means,cov);
            gmmhmm_models=[gmmhmm_models,est_a];
        end
    end
    dist_matrix_dim=GMM_NUM*SEQ_NUM_FOR_EACH_GMM;
    dist_matrix=zeros(dist_matrix_dim);
    for i=1:dist_matrix_dim
        for j=1:i
            [dist]=gmmhmm_dist_likelihood(gmmhmm_models(i),gmmhmm_models(j),500);
            dist_matrix(i,j)=dist;
            dist_matrix(j,i)=dist;
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix,ground_truth_class,['KL \Delta\sigma= ',num2str(factor)]};
    paras{factor_idx}=para;
end
toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
% paras{1}{3}='KL \Delta\sigma= 0.2';
% paras{2}{3}='KL \Delta\sigma= 0.4';
% paras{3}{3}='KL \Delta\sigma= 0.6';
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_likelihood.mat'],'paras');
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_likelihood.mat'])
% draw_multiple_precrecl(paras,'Changing \Sigma','change_sigma_likelihood','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/')


%% R.D. sampling BADMM
factors = [0.2,0.4,0.6];
paras={};
BADMM_SAMPLES=100;
tic
for factor_idx = 1:length(factors)
    para={};
    factor=factors(factor_idx);
    disp(factor);
    GMM_NUM=5;
    SEQ_NUM_FOR_EACH_GMM=10;
    gmmhmm_models=[];
    for i=1:GMM_NUM
        eg.dim=2;
        eg.mu=[[2+i*factor,2+i*factor];[5+i*factor,5+i*factor]]';
        eg.trans_mat=[.8 .2; .2 .8;];
        eg.weights=station_dist(eg.trans_mat);
        r_num=rand(2);
        r_num=(r_num+r_num')/2;
        temp = 0.2*expm(factor*i*r_num);
        eg.covariance=cat(3,temp,temp);
        gmmhmm_eg1 = gmmhmm(eg.dim,eg.trans_mat,eg.weights,eg.mu,eg.covariance);
        n=100; % length of chain
        for j=1:SEQ_NUM_FOR_EACH_GMM
            rng(i*j); % make result reproducible
            data = gmmhmm_eg1.gen_mk_chain(n);
            numst = size(eg.mu,2);
            [a_00, a, gauss]=hmm_fit_jia(data, eg.dim, 1, n,numst,0);
            means=vertcat(gauss.means)';
            cov=cat(3,gauss.sigma);
            est_a = gmmhmm(eg.dim,a,a_00,means,cov);
            gmmhmm_models=[gmmhmm_models,est_a];
        end
    end
    dist_matrix_dim=GMM_NUM*SEQ_NUM_FOR_EACH_GMM;
    dist_matrix_D1=zeros(dist_matrix_dim);
    dist_matrix_D23=zeros(dist_matrix_dim);
    for i=1:dist_matrix_dim
        for j=1:i
            [dist, matching] = gmm_wass_dist_mc_BADMM_posterior(gmmhmm_models(i),gmmhmm_models(j),BADMM_SAMPLES);
            dist_matrix_D1(i,j)=dist;
            dist_matrix_D1(j,i)=dist;
            if any(isnan(matching(:)))
                error('%f',matching);
            end
            D = pdist2(gmmhmm_models(i),gmmhmm_models(j),'gaussianMixture');
            [dist2]=gmmhmm_dist_naive_transmatdist(gmmhmm_models(i),gmmhmm_models(j),matching, D);
            dist_matrix_D23(i,j)=dist2;
            dist_matrix_D23(j,i)=dist2;
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix_D1,dist_matrix_D23,ground_truth_class,['Improved MAW \Delta\sigma= ',num2str(factor)]};
    paras{factor_idx}=para;
end
toc

if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
% paras{1}{4}='Improved MAW \Delta\sigma= 0.2';
% paras{2}{4}='Improved MAW \Delta\sigma= 0.4';
% paras{3}{4}='Improved MAW \Delta\sigma= 0.6';
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_BADMM.mat'],'paras');
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_BADMM.mat'])
% draw_multiple_precrecl(paras,'Changing \Sigma','change_sigma_BADMM','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/',0.1)


%% Compare naive and KL and BADMM
% concatenated_para={};
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_likelihood.mat']);
% concatenated_para{1}=paras;
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_naive.mat']);
% concatenated_para{2}=paras;
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_BADMM.mat']);
% concatenated_para{3}=paras;
% draw_multiple_precrecl_multimethods(concatenated_para,'Changing \sigma','change_sigma_naive_vs_likelihood_vs_BADMM','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/',0.1);

%% Compare naive and KL for 9 figure way
concatenated_para={};
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_likelihood.mat']);
concatenated_para{1}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_naive.mat']);
concatenated_para{2}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_BADMM.mat']);
concatenated_para{3}=paras;
draw_multiple_precrecl({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying \Sigma, \Delta \sigma=0.2','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_2','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/9fig/',0)
draw_multiple_precrecl({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying \Sigma, \Delta \sigma=0.4','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_4','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/9fig/',0)
draw_multiple_precrecl({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying \Sigma, \Delta \sigma=0.6','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_6','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/9fig/',0)

%% Distance of seqs from HMM1 to HMM1 - HMM5
draw_varianceplot({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying \Sigma, \Delta \sigma=0.2','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_2',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying \Sigma, \Delta \sigma=0.4','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_4',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying \Sigma, \Delta \sigma=0.6','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_6',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)

%% Distance of seqs from HMM3 to HMM1 - HMM5
draw_varianceplot_from3({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying \Sigma, \Delta \sigma=0.2','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_2_from3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot_from3({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying \Sigma, \Delta \sigma=0.4','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_4_from3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot_from3({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying \Sigma, \Delta \sigma=0.6','change_sigma_naive_vs_likelihood_vs_BADMM_deltasigma0_6_from3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
