classdef gmm
    %Class of Gaussian Mixture Model
    %   dim: dimension for each gaussian peak;
    %   mode_num: Amount of Gaussian mixture peak in the gmm;
    %   mu: a list, each item is the mu vector (1*dim) for each Gaussian peak;
    %   covariance: a list, each item is the covariance matrix (dim*dim) for each
    %   Gaussian peak;
    % Authors: Yukun Chen, College of Information Science and Technology, Penn State University
    % Contact: cykustc@gmail.com
    
    properties
        dim
        mode_num
        mu
        %mu is stored as:
        %[mu_1_x,mu_2_x,mu_3_x,...]
        %[mu_1_y,mu_2_y,mu_3_y,...]
        %[mu_1_..,mu_2_..,mu_3_..,...]
        covariance %3d cat in the 3rd dimension
        inv_covariance %3d cat in the 3rd dimension
        %covariance is stored as:
        %[cov_1_xx,cov_1_xy,mu_1_xz,...]  .      .
        %[cov_1_yx,cov_1_yy,mu_3_yz,...].     .
        %[cov_1_..x,cov_1_..y,mu_3_..z,...].
        weights % Sum of weights should be 1
    end
    
    methods
        function obj = gmm(dim,mode_num,mu,covariance,weights)
            obj.dim=dim;
            obj.mode_num=mode_num;
            obj.mu=mu;
            obj.covariance=covariance;
            for i=1:size(obj.covariance,3)
               obj.covariance(:,:,i)=make_symmetric(obj.covariance(:,:,i));
            end
            %TRICK: In case covariance(:,:,i) is singlar, add an appropriate diagnol matrix
            for i=1:obj.mode_num
                if(cond(obj.covariance(:,:,i))>10^4)
                    max_eigenvalue=max(eig(obj.covariance(:,:,i)));
                    obj.covariance=obj.covariance+repmat(max_eigenvalue/10^5*eye(obj.dim), 1, 1, 3);
                end
            end
            obj.inv_covariance=zeros(size(covariance));
            for i=1:size(covariance,3)
                obj.inv_covariance(:,:,i)=inv(covariance(:,:,i));
            end
            if (sum(weights)-1 >= 1e-3)
                error('Sum of gmm weights should be 1!');
            end
            obj.weights=weights;
        end
        
        function rndnum = rnd(obj)
            %Gererate 1 rnd number from GMM model
            temprnd=unifrnd(0,1);
            cum_val = 0;
            for i=1:length(obj.weights)
                if temprnd>=cum_val && temprnd<=cum_val+obj.weights(i)
                    rndnum=mvnrnd(double(obj.mu(:,i)),obj.covariance(:,:,i));
                    break
                else
                    cum_val = cum_val+obj.weights(i);
                end
            end
        end
        
        function [rndnum, component] = rnd_component_pair_gen(obj)
            %Gererate 1 rnd number from GMM model
            temprnd=unifrnd(0,1);
            cum_val = 0;
            for i=1:length(obj.weights)
                if temprnd>=cum_val & temprnd<=cum_val+obj.weights(i)
                    rndnum=mvnrnd(double(obj.mu(:,i)),obj.covariance(:,:,i));
                    component = i;
                    break
                else
                    cum_val = cum_val+obj.weights(i);
                end
            end
        end
        
        function rndnumlist = rndlist(obj, n)
            % Generate a vector of n items, each is a random vector of
            % obj.dim dimmension from GMM model.
            sample_num_in_each_mode = mnrnd(n,obj.weights); % sample number in each mode
            rndnumlist=zeros(obj.dim,n); % vector to store the final results, initialized with 0s.
            %         [r1_1, r2_1,...,rn_1]
            %obj.dim  [r1_1, r2_1,...,rn_1]
            %         [r1_1, r2_1,...,rn_1]
            %                   n
            initial_index=1;
            for i = 1:length(sample_num_in_each_mode)
%                 disp(i)
%                 disp(obj.mu(:,i));
%                 disp(obj.covariance(:,:,i));
                i_th_sample_size = sample_num_in_each_mode(i);
                temp_cov = (obj.covariance(:,:,i) + obj.covariance(:,:,i)')/2;
                rndnums_in_ith_mode=mvnrnd(double(obj.mu(:,i)),temp_cov,i_th_sample_size);
                if i_th_sample_size>0
%                     disp(size(rndnums_in_ith_mode'))
%                     disp(size(rndnumlist(:,initial_index:initial_index+i_th_sample_size-1)))
                    rndnumlist(:,initial_index:initial_index+i_th_sample_size-1)=rndnums_in_ith_mode';
                    initial_index = initial_index + i_th_sample_size;
                end
            end
        end
        function [rndnumlist, complist] = rndlist_complist_gen(obj, n)
            % Generate 2 vectors
            % ist vector is of n items, each is a random vector of
            % obj.dim dimmension from GMM model.
            % 2nd vector is of n itmes, each denotes which components it
            % belongs to.
            sample_num_in_each_mode = mnrnd(n,obj.weights); % sample number in each mode
            rndnumlist=zeros(obj.dim,n); % vector to store the final results, initialized with 0s.
            %         [r1_1, r2_1,...,rn_1]
            %obj.dim  [r1_1, r2_1,...,rn_1]
            %         [r1_1, r2_1,...,rn_1]
            %                   n
            complist=zeros(1,n);
            initial_index=1;
            for i = 1:length(sample_num_in_each_mode)
%                 disp(i)
%                 disp(obj.mu(:,i));
%                 disp(obj.covariance(:,:,i));
                i_th_sample_size = sample_num_in_each_mode(i);
                temp_cov = (obj.covariance(:,:,i) + obj.covariance(:,:,i)')/2;
                rndnums_in_ith_mode=mvnrnd(double(obj.mu(:,i)'),temp_cov,i_th_sample_size);
                if i_th_sample_size>0
%                     disp(size(rndnums_in_ith_mode'))
%                     disp(size(rndnumlist(:,initial_index:initial_index+i_th_sample_size-1)))
                    rndnumlist(:,initial_index:initial_index+i_th_sample_size-1)=rndnums_in_ith_mode';
                    complist(initial_index:initial_index+i_th_sample_size-1)=i;
                    initial_index = initial_index + i_th_sample_size;
                end
            end
        end
        function [rndnumlist, complist] = rndlist_complist_gen_depreciated(obj, n)
            % Generate 2 vectors rndnumlist and complist,
            % rndnumlist is a vector of n items, each is a random vector of
            % obj.dim dimmension from GMM model.
            % complist is a vector each record the component, from which the corresponding
            % item in rndnumlist is generated from.
            rndnumlist=zeros(obj.dim,n);
            complist=zeros(1,n);
            for i=1:n
                [rndnumlist(:,i),complist(i)] = obj.rnd_component_pair_gen(); 
            end
        end
        
        function plot_rndsamples(obj,sample_num)
            if obj.dim==1
                % hist(obj.rndlist(sample_num),binnum);
                [rndnumlist, complist] = obj.rndlist_complist_gen(sample_num);
                scatter(rndnumlist,zeros(1,length(rndnumlist)),50,complist)
            elseif obj.dim==2
                [rndnumlist, complist] = obj.rndlist_complist_gen(sample_num);
                scatter(rndnumlist(1,:),rndnumlist(2,:),50,complist)
            elseif obj.dim==3
                [rndnumlist, complist] = obj.rndlist_complist_gen(sample_num);
                scatter3(rndnumlist(1,:),rndnumlist(2,:),rndnumlist(3,:),50,complist)
            else
                disp('rnd generator for dim>1 has not been implemented!')
            end
        end
        
        %%Depreciated:
        function rndnumlist = rndlist2(obj, n)
            % Generate a vector of n items, each is a random vector of
            % obj.dim dimmension from GMM model.
            rndnumlist=zeros(obj.dim,n);
            for i=1:n
                rndnumlist(:,i) = obj.rnd();
            end
        end
        
        function softmembership = softmembership(obj, points)
            n=size(points,2);
            p=gmdistribution(obj.mu',obj.covariance,obj.weights);
            softmembership=p.posterior(points');
            softmembership=softmembership';
%             softmembership = zeros(obj.mode_num, n);
%             for i=1:n
%                 tmp=p.posterior(points(:,i)');
%                 softmembership(:,i)=tmp';
%             end
        end
    end
end

