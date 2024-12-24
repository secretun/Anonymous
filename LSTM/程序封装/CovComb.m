function Cov = CovComb(X,step,overlap)
    [row,col,num] = size(X);
    cov_count = 0;
    n = floor((col-step)/(step-overlap))+1;
    Cov = zeros(row,row,num*n);
    for i = 1:num
        for j = 1:(step-overlap):col
            if j+step-1 <= col
                cov_count = cov_count + 1;  
                seg = X(:,j:(j+step-1),i);
                Cov(:,:,cov_count) = covariances(seg);
            end
        end
    end
end