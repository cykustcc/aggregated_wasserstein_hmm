classdef gmmhmm < gmm
    % Class of GMM-HMM Model

    properties
        trans_mat
    end
    
    methods
        function obj=gmmhmm(dim,trans_mat,station_dist,mu,covariance)
            mode_num=length(trans_mat);
            weights=station_dist;            
            obj@gmm(dim,mode_num,mu,covariance,weights);
            obj.trans_mat=trans_mat;
        end
        function [chain,seq] = gen_mk_chain(obj,n)
            chain = zeros(obj.dim,n);
            EMIS = eye(length(obj.trans_mat));
            [seq,~] = hmmgenerate(n,obj.trans_mat,EMIS);
            for i = 1:obj.mode_num
                chain(:,seq==i)=mvnrnd(obj.mu(:,i)',obj.covariance(:,:,i),length(seq(seq==i)))';
            end
        end
%         function gmm_model = gmm(obj)
%             gmm_model = gmm(obj.dim,obj.mode_num,obj.mu,obj.covariance,obj.weights);
%         end
        function print(obj)
            disp('dim:');
            disp(obj.dim);
            disp('weights:');
            disp(obj.weights);
            disp('tran_mat:');
            disp(obj.trans_mat);
            disp('means');
            disp(obj.mu);
            disp('covariances:');
            disp(obj.covariance);
        end
        
        function obj=get_rid_of_weight_0_component(obj)
            % If the weight of one component is 0, get rid of that
            % component and corresponding row/column in trans_mat
            
            index=1:obj.mode_num;
            index(abs(obj.weights(:)-0)<10^-5)=[];

            if length(index)<obj.mode_num
                obj.mode_num=size(index,2);
                obj.weights=obj.weights(index);
                obj.mu=obj.mu(:,index);
                obj.trans_mat=obj.trans_mat(index,index);
                obj.covariance=obj.covariance(:,:,index);
            end
        end
    end
    
end

