function proj = projection3D( A, B, C, P, mode)

%{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Local coordinate extraction of P, on a coordinate system defined by A, B
and C

A, B, C and P are [N,3] matrices (N = number of frames)

optional parameter:
- mode: if mode = 0 => projection; if mode = 1 => inverse projection

%}

if nargin < 5
    mode=0;
end

proj = ones(size(P))*nan;

nans=isnan(A) | isnan(B) | isnan(C);

%Coordinate system (vectorized version)
v1 = B - A;
v1 = normr(v1);
v2 = cross(v1, C - A, 2);
v2 = normr(v2);
v3 = cross(v2, v1, 2)';
v2 = v2';
v1 = v1';

%(vectorized for speed)
Mat(1,:,:) = v1;
Mat(2,:,:) = v2;
Mat(3,:,:) = v3;
Mat = permute(Mat,[2 1 3]);
PminusA = P-A;

%Projection at each time frame
for i=1:size(P,1)
        
    %faster
    projmat = Mat(:,:,i);
    
    if(mode == 0) %normal projection
        proj(i,:) = (PminusA(i,:))*projmat;
    end
    if(mode == 1) %inverse projection
        if ~isnan(sum(projmat(:))) %no nan data
            %invprojmat = pinv(projmat);
            %proj(i,:) = P(i,:)*invprojmat + A(i,:);
            %faster (if compatible with your matlab version)
            if(abs(det(projmat)) > eps) %check if matrix is nearly singular
                proj(i,:) = P(i,:)/projmat + A(i,:);
            else %slower but better if matrix is nearly singular
                proj(i,:) = P(i,:)*pinv(projmat) + A(i,:);
            end
        end
    end
end
proj(nans)=nan;

end