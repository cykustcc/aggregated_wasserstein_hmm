%%%%%%%%%%%%%%%%%%%
%% Perturb trans_mat
%%%%%%%%%%%%%%%%%%%
%% Test trans_mat generator
factors = [0.3,0.6,0.9];
sigma=[1,0;0,1]
base_trans_mat=[40,10;10,40];

%% gmm_hmm_naive
paras={};

rng(2)
trans_mat_bank={};
for i=1:GMM_NUM
    trans_mat_bank{i}=gen_trans_mat(base_trans_mat,0.1)
end
tic
base_trans_matrix=[0.8,0.2;0.2,0.8];
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
        factor = factors(factor_idx);
%         a = expm(factor*rand(2));
%         eg.trans_mat = bsxfun(@times,a,1./sum(a,2));
        
        eg.trans_mat=(1-factor)*base_trans_matrix+factor*trans_mat_bank{i};
        eg.weights=station_dist(eg.trans_mat);
        eg.covariance=cat(3,sigma,sigma);
        gmmhmm_eg1 = gmmhmm(eg.dim,eg.trans_mat,eg.weights,eg.mu,eg.covariance);
        n=1000; % length of chain
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
            [dist,matching]=gmm_wass_dist_naive(gmmhmm_models(i),gmmhmm_models(j), D);
            dist2 = gmmhmm_dist_naive_transmatdist(gmmhmm_models(i),gmmhmm_models(j),matching, D);
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
    para={dist_matrix,dist_matrix2,ground_truth_class,['MAW \Delta t= ',num2str(factor)]};
    paras{factor_idx}=para;
end
toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
% paras{1}{4}='MAW \Delta t= 0.3';
% paras{2}{4}='MAW \Delta t= 0.6';
% paras{3}{4}='MAW \Delta t= 0.9';
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_naive.mat'],'paras');
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_naive.mat'])
% draw_multiple_precrecl(paras,'Changing Trans. Mat. ','change_transmat_naive','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/')


%% KL divergence likelihood --Juang's paper
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
        factor = factors(factor_idx);
%         a = expm(factor*rand(2));
%         eg.trans_mat = bsxfun(@times,a,1./sum(a,2));
        eg.trans_mat=(1-factor)*base_trans_matrix+factor*trans_mat_bank{i};
        eg.weights=station_dist(eg.trans_mat);
        eg.covariance=cat(3,sigma,sigma);
        gmmhmm_eg1 = gmmhmm(eg.dim,eg.trans_mat,eg.weights,eg.mu,eg.covariance);
        n=1000; % length of chain
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
    para={dist_matrix,ground_truth_class,['KL \Delta t= ',num2str(factor)]};
    paras{factor_idx}=para;
end
toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
% paras{1}{3}='KL \Delta t= 0.3';
% paras{2}{3}='KL \Delta t= 0.6';
% paras{3}{3}='KL \Delta t= 0.9';
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_likelihood.mat'],'paras');
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_likelihood.mat'])
% draw_multiple_precrecl(paras,'Changing Trans. Mat ','change_transmat_likelihood','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/')


%% Sampling BADMM registration distance.
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
        eg.mu=[[2,2];[5,5]]';
        factor = factors(factor_idx);
%         a = expm(factor*rand(2));
%         eg.trans_mat = bsxfun(@times,a,1./sum(a,2));
        eg.trans_mat=(1-factor)*base_trans_matrix+factor*trans_mat_bank{i};
        eg.weights=station_dist(eg.trans_mat);
        eg.covariance=cat(3,sigma,sigma);
        gmmhmm_eg1 = gmmhmm(eg.dim,eg.trans_mat,eg.weights,eg.mu,eg.covariance);
        n=1000; % length of chain
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
    para={dist_matrix_D1,dist_matrix_D23,ground_truth_class,['Improved MAW \Delta t= ',num2str(factor)]};
    paras{factor_idx}=para;
end
toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
% paras{1}{4}='Improved MAW \Delta t= 0.3';
% paras{2}{4}='Improved MAW \Delta t= 0.6';
% paras{3}{4}='Improved MAW \Delta t= 0.9';
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_BADMM.mat'],'paras');
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_BADMM.mat'])
% draw_multiple_precrecl(paras,'Changing T','change_transmat_BADMM','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/',0.2)


%% Compare naive and KL and BADMM
% concatenated_para={};
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_likelihood.mat']);
% concatenated_para{1}=paras;
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_naive.mat']);
% concatenated_para{2}=paras;
% load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_BADMM.mat']);
% concatenated_para{3}=paras;
% draw_multiple_precrecl_multimethods(concatenated_para,'Changing Trans. Mat.','change_transmat_naive_vs_likelihood_vs_BADMM','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/',1);

%% Compare naive and KL for 9 figure way
concatenated_para={};
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_likelihood.mat']);
concatenated_para{1}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_naive.mat']);
concatenated_para{2}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_BADMM.mat']);
concatenated_para{3}=paras;
draw_multiple_precrecl({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying Trans. Mat., \Delta t=0.3','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_3','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/9fig/',1)
draw_multiple_precrecl({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying Trans. Mat., \Delta t=0.6','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_6','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/9fig/',1)
draw_multiple_precrecl({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying Trans. Mat., \Delta t=0.9','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_9','/Users/MAC/Dropbox/GMMHMM/src/','imgs/perturbaion_exp/9fig/',1)

draw_varianceplot({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying Trans. Mat., \Delta t=0.3','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_3','/Users/MAC/Dropbox/GMMHMM/src/','imgs/varianceplot/9fig/',1)
draw_varianceplot({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying Trans. Mat., \Delta t=0.6','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_6','/Users/MAC/Dropbox/GMMHMM/src/','imgs/varianceplot/9fig/',1)
draw_varianceplot({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying Trans. Mat., \Delta t=0.9','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_9','/Users/MAC/Dropbox/GMMHMM/src/','imgs/varianceplot/9fig/',1)

