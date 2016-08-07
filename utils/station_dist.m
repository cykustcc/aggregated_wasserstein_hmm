function [st,D] = station_dist( TRANS )
    % return the stationary distribution [1*length(TRANS)]
    % of the Markov Chain defined by the transition matrix TRANS.
    [V,D] = eig(TRANS'); %each col of V: eigenvector, D(i,i) is eigenvalue.
    st=V(:,1);
    st=st'/sum(st);
end

