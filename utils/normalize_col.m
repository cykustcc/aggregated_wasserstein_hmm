function matching_c=normalize_col(matching,gmmhmm2_weights)
    matching_c=bsxfun(@times,matching,1./gmmhmm2_weights);
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