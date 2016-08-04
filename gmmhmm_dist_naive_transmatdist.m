function [Wassmat_dist] = gmmhmm_dist_naive_transmatdist(gmmhmm1,gmmhmm2,matching, D)
%GMMHMM_DIST_NAIVE naive way of computing gmmhmm distance
    matching=get_rid_of_weight_all0rowcol_matching(matching,gmmhmm1,gmmhmm2);
    gmmhmm1=gmmhmm1.get_rid_of_weight_0_component();
    gmmhmm2=gmmhmm2.get_rid_of_weight_0_component();
%     disp(matching);
%     disp(gmmhmm1.weights);
%     disp(gmmhmm2.weights);
    matching_c=normalize_col(matching, gmmhmm2.weights);
    matching_r=normalize_row(matching, gmmhmm1.weights);
    trans_mat_1to2=matching_c'*gmmhmm1.trans_mat*matching_r;
    trans_mat_2to1=matching_r*gmmhmm2.trans_mat*matching_c';
    Wassmat_dist=0;
    for i=1:gmmhmm2.mode_num
        if abs(sum(trans_mat_1to2(i,:))-1)>10^-5            
%             disp('problem row:')
%             disp(trans_mat_1to2(i,:));
%             disp(sum(trans_mat_1to2(i,:)));
%             disp('weights:');
%             disp(gmmhmm1.weights);
%             disp(sum(gmmhmm1.weights));
%             disp(gmmhmm2.weights);
%             disp(sum(gmmhmm2.weights));
%             disp('matching:');
%             disp(matching);
%             disp(sum(sum(matching)));
%             disp(matching_c);
%             disp(sum(matching_c,1));
%             disp(matching_r);
%             disp(sum(matching_r,2));
%             disp('old trans mat:');
%             disp(gmmhmm1.trans_mat);
%             disp(sum(gmmhmm1.trans_mat,1));
%             disp(sum(gmmhmm1.trans_mat,2));
%             disp(gmmhmm2.trans_mat);
%             disp(sum(gmmhmm2.trans_mat,1));
%             disp(sum(gmmhmm2.trans_mat,2));
%             disp('new trans mat:');
%             
%             disp(trans_mat_1to2);
%             disp(sum(trans_mat_1to2,1));
%             disp(sum(trans_mat_1to2,2));
%             disp(trans_mat_2to1);
%             disp(sum(trans_mat_2to1,1));
%             disp(sum(trans_mat_2to1,2));            
% 	    disp('covariance mat:')
% 	    disp(gmmhmm1.covariance);
% 	    disp(gmmhmm2.covariance);
%         
        trans_mat_1to2(i,:)=trans_mat_1to2(i,:)/sum(trans_mat_1to2(i,:));
        end
        
        if nargin < 4
            gmm_trans_mat_1_i=gmm(gmmhmm2.dim,gmmhmm2.mode_num,gmmhmm2.mu,gmmhmm2.covariance,trans_mat_1to2(i,:));
            gmm_trans_mat_2_i=gmm(gmmhmm2.dim,gmmhmm2.mode_num,gmmhmm2.mu,gmmhmm2.covariance,gmmhmm2.trans_mat(i,:));
            [mat_dist, x, lambda] = gmm_wass_dist_naive(gmm_trans_mat_1_i,gmm_trans_mat_2_i);
        else
            [~, mat_dist, x, lambda] = OptimalTransport_LP(D, trans_mat_1to2(i,:)',gmmhmm2.trans_mat(i,:)');
        end
        Wassmat_dist=Wassmat_dist+mat_dist*gmmhmm2.weights(i);
    end
    for i=1:gmmhmm1.mode_num
        if abs(sum(trans_mat_2to1(i,:))-1)>10^-5
%             disp(trans_mat_2to1(i,:));
%             disp(sum(trans_mat_1to2(i,:)));
%             disp('weights:');
%             disp(gmmhmm1.weights);
%             disp(sum(gmmhmm1.weights));
%             disp(gmmhmm2.weights);
%             disp(sum(gmmhmm2.weights));
%             disp('matching:');
%             disp(matching);
%             disp(sum(sum(matching)));
%             disp(matching_c);
%             disp(sum(matching_c,1));
%             disp(matching_r);
%             disp(sum(matching_r,2));
%             disp('old trans mat:');
%             disp(gmmhmm1.trans_mat);
%             disp(sum(gmmhmm1.trans_mat,1));
%             disp(sum(gmmhmm1.trans_mat,2));
%             disp(gmmhmm2.trans_mat);
%             disp(sum(gmmhmm2.trans_mat,1));
%             disp(sum(gmmhmm2.trans_mat,2));
%             disp('new trans mat:');
%             disp(trans_mat_1to2);
%             disp(sum(trans_mat_1to2,1));
%             disp(sum(trans_mat_1to2,2));
%             disp(trans_mat_2to1);
%             disp(sum(trans_mat_2to1,1));
%             disp(sum(trans_mat_2to1,2));
%             
%             disp('covariance mat:')
%             disp(gmmhmm1.covariance);
%             disp(gmmhmm2.covariance);
            
            trans_mat_2to1(i,:)=trans_mat_2to1(i,:)/sum(trans_mat_2to1(i,:));
        end
        if nargin < 4
            gmm_trans_mat_1_i=gmm(gmmhmm1.dim,gmmhmm1.mode_num,gmmhmm1.mu,gmmhmm1.covariance,gmmhmm1.trans_mat(i,:));
            gmm_trans_mat_2_i=gmm(gmmhmm1.dim,gmmhmm1.mode_num,gmmhmm1.mu,gmmhmm1.covariance,trans_mat_2to1(i,:));
        
            [mat_dist, x, lambda] = gmm_wass_dist_naive(gmm_trans_mat_1_i,gmm_trans_mat_2_i);
        else
            [~, mat_dist, x, lambda] = OptimalTransport_LP(D,gmmhmm1.trans_mat(i,:)',trans_mat_2to1(i,:)');
        end
        Wassmat_dist=Wassmat_dist+mat_dist*gmmhmm1.weights(i);
    end
end

function matching=get_rid_of_weight_all0rowcol_matching(matching,gmmhmm1,gmmhmm2)
    % If the sum of the row or col is 0, get rid of that row/col
    index1=[];
    index2=[];
    flag=0;
    for i=1:gmmhmm1.mode_num
        if abs(gmmhmm1.weights(i)-0)>10^-5
            index1=[index1,i];
        else
            flag=1;
        end
    end
    for i=1:gmmhmm2.mode_num
        if abs(gmmhmm2.weights(i)-0)>10^-5
            index2=[index2,i];
        else
            flag=1;
        end
    end
    if flag==1
        matching=matching(index1,index2);
    end
end
        
% function matching_c=normalize_col(matching,gmmhmm2)
%     matching_c=bsxfun(@times,matching,1./gmmhmm2.weights);
% %     matching_c=zeros(size(matching));
% %     for i=1:m
% %         if gmmhmm2.weights(i)~=0
% %             matching_c(:,i)=matching(:,i)/gmmhmm2.weights(i);
% %         end
% %     end
% end
% 
% function matching_r=normalize_row(matching,gmmhmm1)
%     matching_r=bsxfun(@times,matching,1./gmmhmm1.weights');
% %     matching_r=zeros(size(matching));
% %     for i=1:n
% %         if gmmhmm1.weights(i)~=0
% %             matching_r(i,:)=matching(i,:)/gmmhmm1.weights(i);
% %         end
% %     end
% end
