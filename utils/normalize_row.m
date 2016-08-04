function matching_r=normalize_row(matching,gmmhmm1_weights)
    matching_r=bsxfun(@times,matching,1./gmmhmm1_weights');
end
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