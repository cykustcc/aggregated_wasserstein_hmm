%%%%%%%%%%%%%%%%%%%
%% Perturb sigma
%%%%%%%%%%%%%%%%%%%
%% gmm_hmm_naive
factors = [0.2,0.4,0.6];
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
            dist = MAW(gmmhmm_models(i),gmmhmm_models(j), 1);
            dist_matrix(i,j)=dist;
            dist_matrix(j,i)=dist; 
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix,ground_truth_class,['MAW \Delta\sigma= ',num2str(factor)]};
    paras{factor_idx}=para;
end
% toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end

save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_MAW.mat'],'paras');


%% KL divergence likelihood --Juang's paper
factors = [0.2,0.4,0.6];
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
            [dist]=gmmhmm_KL(gmmhmm_models(i),gmmhmm_models(j),500);
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
% toc
if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_KL.mat'],'paras');


%% R.D. sampling BADMM
factors = [0.2,0.4,0.6];
paras={};
SAMPLE_SIZE=100;
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
            options.method = 'sinkhorn';
            options.sample_size = SAMPLE_SIZE;
            dist = IAW(gmmhmm_models(i),gmmhmm_models(j), 1, options);
            dist_matrix(i,j)=dist;
            dist_matrix(j,i)=dist;
        end
    end
    ground_truth_class=zeros(1,dist_matrix_dim);
    for k=1:GMM_NUM
      ground_truth_class(1+(k-1)*SEQ_NUM_FOR_EACH_GMM:k*SEQ_NUM_FOR_EACH_GMM) = k;
    end
    para={dist_matrix,ground_truth_class,['IAW \Delta\sigma= ',num2str(factor)]};
    paras{factor_idx}=para;
end
% toc

if ~exist([gmmhmm_projectroot,'/data/perturbaion_exp/'], 'dir')
  mkdir([gmmhmm_projectroot,'/data/perturbaion_exp/']);
end
save([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_IAW.mat'],'paras');


%% Load stored MAW, KL and IAW distances.
concatenated_para={};
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_KL.mat']);
concatenated_para{1}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_MAW.mat']);
concatenated_para{2}=paras;
load([gmmhmm_projectroot,'/data/perturbaion_exp/','change_sigma_IAW.mat']);
concatenated_para{3}=paras;

%% Distance of seqs from HMM1 to HMM1 - HMM5
draw_varianceplot({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying \Sigma, \Delta \sigma=0.2','change_sigma_MAW_vs_KL_vs_IAW_deltasigma0_2',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying \Sigma, \Delta \sigma=0.4','change_sigma_MAW_vs_KL_vs_IAW_deltasigma0_4',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying \Sigma, \Delta \sigma=0.6','change_sigma_MAW_vs_KL_vs_IAW_deltasigma0_6',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)

%% Distance of seqs from HMM3 to HMM1 - HMM5
draw_varianceplot_from3({concatenated_para{1}{1},concatenated_para{2}{1},concatenated_para{3}{1}},'Varying \Sigma, \Delta \sigma=0.2','change_sigma_MAW_vs_KL_vs_IAW_deltasigma0_2_from3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot_from3({concatenated_para{1}{2},concatenated_para{2}{2},concatenated_para{3}{2}},'Varying \Sigma, \Delta \sigma=0.4','change_sigma_MAW_vs_KL_vs_IAW_deltasigma0_4_from3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
draw_varianceplot_from3({concatenated_para{1}{3},concatenated_para{2}{3},concatenated_para{3}{3}},'Varying \Sigma, \Delta \sigma=0.6','change_sigma_MAW_vs_KL_vs_IAW_deltasigma0_6_from3',gmmhmm_projectroot,'/imgs/varianceplot/9fig/',0)
