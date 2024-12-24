function result = threeclassconfusion(YTest,YPred)
    if iscategorical(YTest)
        YTest1 = string(YTest);
        YTest1 = double(YTest1);
    end
    
    if iscategorical(YPred)
        YPred1 = string(YPred);
        YPred1 = double(YPred1);
    end
    
    labels = unique(YTest1);
    n = length(labels);
    result = zeros(n,n);
    for i = 1:n
        for j = 1:n
            idx = (YTest1==labels(i));
            YPred2 = YPred1(idx);
            result(i,j) = sum(YPred2==labels(j))/length(YPred2);
        end
    end
        
    
end