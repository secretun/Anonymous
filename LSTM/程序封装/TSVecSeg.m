function result = TSVecSeg(X,step) 
% Tangent Space Vector Segment
    vec_count = 0;
    [~,col] = size(X);
    for i = 1:step:col
        vec_count = vec_count + 1;
        A = real(X(:,i:(i+step-1)));
        [r1,c1] = size(A);
        B = reshape(A,1,r1*c1);
        B = B - mean(B);
        B = B ./ std(B);
        A = reshape(B,r1,c1);
        result{vec_count,1} = A;
    end
end