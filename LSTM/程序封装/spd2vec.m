function Feat = spd2vec(X)
    % X is a symmetric positive matrix
    N_elec = size(X,1);
    NTrial = size(X,3);
    index = reshape(triu(ones(N_elec)),N_elec*N_elec,1)==1;
    for i=1:NTrial
        Tn = X(:,:,i);
        tmp = reshape(sqrt(2)*triu(Tn,1)+diag(diag(Tn)),N_elec*N_elec,1);
        Feat(:,i) = tmp(index);   
    end
end