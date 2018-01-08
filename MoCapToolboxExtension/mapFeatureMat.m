function out = mapFeatureMat(X, degree)

%{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Multinomial products of (1 to degree) features of X with repetition
%}


if nargin < 2
    degree = 2;
end


dim = size(X,2);
out = [];
for i=1:degree
    combs = nmultichoosek(1:dim,i);
    for c = 1:size(combs,1)
        out(:,end+1) = prod(X(:,combs(c,:)),2);
    end
    
end

end

function combs = nmultichoosek(values, k)
%(combination with repetition)
%// Return number of multisubsets or actual multisubsets.
if numel(values)==1
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
end