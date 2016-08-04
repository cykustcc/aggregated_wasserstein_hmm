function [precT, reclT, prec, recl, cutof, accu] = precisionRecall( score, label, varargin )
% precisionRecallPlot plot the precision-recall curve
% Input:
%   score -- the ranking measure
%   label -- friendship label 0/1
%   varargin -- other arguments controlling the plot
% Output:
%   precT -- precision values on given query 
%   reclT -- given values of recall as the query
%   line  -- the line object of plot
%   prec  -- the precision used for plotting
%   recl  -- the recall used for plotting
%   cutof -- the cutoff score for the recall-precision calculation
%   accu  -- the accuracy

    [prec, recl, cutof, accu] = precisionRecall( score, label );
    
    precT = zeros(3,1);
    reclT = [0.3, 0.5, 0.7];
    for j = 1:length(reclT)
        for tt = 1:length(recl)
            if recl(tt) == reclT(j)
                precT(j) = prec(tt);
            else if recl(tt) < reclT(j) && recl(tt+1) > reclT(j)
                    precT(j) = interpolate(recl(tt), recl(tt+1), prec(tt), prec(tt+1), reclT(j));
                    break;
                end
            end
        end
    end


    function [prec, recl, cutoff, accu] = precisionRecall( score, label )

        if length(score) ~= length(label)
            error('length of score and label does not match.')
        end

        n = length(score);
        data = zeros(n, 2);
        data(:,1) = score;
        data(:,2) = label;
                
        % randomize the vector first
        ind = randperm( length(score) );
        data = data(ind,:);

        
        [~, ind] = sort(data(:,1), 'descend');
        data = data(ind, :);

        step = max( round(n / 100), 1);
        totalPos = sum(data(:,2)==1);

        prec = zeros(1,1);
        recl = zeros(1,1);
        cutoff = zeros(1,1);
        accu = zeros(1,1);

        ind = 1;
        for i = 1:step:n
            head = data(1:i,:);
            tail = data(i+1:end,:);
            npos = sum(head(:,2)==1);   % number of true positive
            nneg = sum(tail(:,2)==0);   % number of true negative
            prec(ind) = npos / i;
            recl(ind) = npos / totalPos;
            cutoff(ind) = head(end,1);
            accu(ind) = (npos + nneg) / n;
            ind = ind + 1;
        end
    end
    
    function [prec] = interpolate( reca, recb, prea, preb, rect )
        s = (rect - reca) / (recb - reca);
        prec = prea + s * (preb - prea);
    end
end