% pcsquash - compress data using Principal Component Analysis (PCA)
%            into a principal component subspace.  To project back 
%            into the original channel space, use pcexpand() or exproj()
% Usage: 
%          >> [eigenvectors,eigenvalues] = pcsquash(data,ncomps);
%    or    >> [eigenvectors,eigenvalues,compressed,datamean] ...
%                                                    = pcsquash(data,ncomps);
% Inputs:
%          data    - (chans,frames) each row is a channel, each column a time point
%          ncomps  - numbers of components to retain
% Outputs: 
%          eigenvectors = square matrix of (column) eigenvectors 
%          eigenvalues  = vector of associated eigenvalues 
%          compressed   = data compressed into space of the ncomps eigenvectors
%                          with largest eigenvalues (ncomps,frames)
%          Note:        >> compressed = eigenvectors(:,1:ncomps)'*data;
%          datamean     = input data channel (row) means (used internally)

% Tzyy-Ping Jung & Scott Makeig CNL / Salk Institute, La Jolla CA 6-97

function [EigenVectors,EigenValues,Compressed,Datamean]=pcsquash(matrix,ncomps)

if nargin < 1 
   help pcsquash
   return
end
if nargin < 2
   ncomps = 0;
end
if ncomps == 0
  ncomps = size(matrix,1);
end
if ncomps < 1
   help pcsquash
   return
end

data = matrix';                    % transpose data
[n,p]=size(data);                  % now p chans,n time points
if ncomps > p
   fprintf('pcsquash(): components must be <= number of data rows (%d).\n',p);
   help pcsquash
end

Datamean = mean(data);  % remove column (channel) means
data = data-ones(n,1)*Datamean;    % remove column (channel) means
out=data'*data/n;
[V,D] = eig(out);                  % get eigenvectors/eigenvalues
diag(D);
[eigenval,index] = sort(diag(D));
index=rot90(rot90(index));
EigenValues=rot90(rot90(eigenval))';
EigenVectors=V(:,index);

if nargout >= 3
   Compressed = EigenVectors(:,1:ncomps)'*matrix;
end
