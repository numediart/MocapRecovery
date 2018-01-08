function X = naninterp(X,method) 
% Interpolate over NaNs 
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), method); 
return 