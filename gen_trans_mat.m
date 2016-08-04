function trans_mat = gen_trans_mat( base_trans_mat, delta_T )
    n = size(base_trans_mat,2);
    trans_mat=zeros(n);
    for i=1:n
        a=base_trans_mat(i,:);
        trans_mat(i,:)=drchrnd(a*delta_T,1);
    end
end

