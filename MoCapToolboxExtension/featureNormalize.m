function [X_norm, mu, sigma] = featureNormalize(X)
%FEATURENORMALIZE Normalizes the features in X 
%   FEATURENORMALIZE(X) returns a normalized version of X where
%   the mean value of each feature is 0 and the standard deviation
%   is 1.

%If X has only one frame
if size(X,1) == 1
   sigma = ones(size(X,2),1);
   mu = X;
   X_norm = X*0;
   
   return;
   
end

mu = nanmean(X);
X_norm = bsxfun(@minus, X, mu);
sigma = nanstd(X_norm);
X_norm = bsxfun(@rdivide, X_norm, sigma);
X_norm(:,sigma==0)=0;


% ============================================================

end
