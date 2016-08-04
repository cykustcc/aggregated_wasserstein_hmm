function [ matrix ] = make_symmetric( matrix )
    ITER=10;
    for i=1:ITER
        matrix=(matrix+matrix')/2;
    end
end

