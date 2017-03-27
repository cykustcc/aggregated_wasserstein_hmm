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
% tic
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
            dist = MAW(gmmhmm_models(i),gmmhmm_models(j), 0);
            dist_matrix(i,j)=dist;
            dist_matrix(j,i)=dist; 
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix,ground_truth_class,['MAW \Delta t= ',num2str(factor)]};
    paras{factor_idx}=para;
end
% toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end

save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_MAW.mat'],'paras');


%% KL divergence likelihood --Juang's paper
paras={};
% tic
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
            [dist]=gmmhmm_KL(gmmhmm_models(i),gmmhmm_models(j),500);
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
% toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end

save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_KL.mat'],'paras');


%% IAW
paras={};
SAMPLE_SIZE=200;
% tic
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
            options.method = 'sinkhorn';
            options.sample_size = SAMPLE_SIZE;
            dist = IAW(gmmhmm_models(i),gmmhmm_models(j), 0, options);
            dist_matrix(i,j)=dist;
            dist_matrix(j,i)=dist;
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix,ground_truth_class,['IAW \Delta t= ',num2str(factor)]};
    paras{factor_idx}=para;
end
% toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end

save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_IAW.mat'],'paras');


%% Compare MAW, KL and IAW for 9 figure way
concatenated_para={};
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_KL.mat']);
concatenated_para{1}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_MAW.mat']);
concatenated_para{2}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_transmat_IAW.mat']);
concatenated_para{3}=paras;
draw_multiple_precrecl({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying Trans. Mat., \Delta t=0.3','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_3',gmmhmm_projectroot,'/imgs/perturbaion_exp/9fig/',1)
draw_multiple_precrecl({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying Trans. Mat., \Delta t=0.6','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_6',gmmhmm_projectroot,'/imgs/perturbaion_exp/9fig/',1)
draw_multiple_precrecl({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying Trans. Mat., \Delta t=0.9','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_9',gmmhmm_projectroot,'/imgs/perturbaion_exp/9fig/',1)

draw_varianceplot({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying Trans. Mat., \Delta t=0.3','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',1)
draw_varianceplot({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying Trans. Mat., \Delta t=0.6','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_6',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',1)
draw_varianceplot({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying Trans. Mat., \Delta t=0.9','change_transmat_naive_vs_likelihood_vs_BADMM_deltat0_9',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',1)

