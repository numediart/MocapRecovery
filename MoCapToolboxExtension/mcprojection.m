function [proj] = mcprojection(m1, m2, m3, p, mode)

%{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Local coordinate extraction of p, on a coordinate system defined by m1, m2
and m3

m1, m2, m3 and p are MoCap structure with one marker (they can be extracted 
with the function mcgetmarker)

optional parameter:
- mode: if mode = 0 => projection; if mode = 1 => inverse projection

%}

if nargin < 5
    mode = 0;
end

A = m1.data;
B = m2.data;
C = m3.data;
proj = p;
proj.data=nan*proj.data;
P = p.data;

nans=isnan(A) | isnan(B) | isnan(C);

%Coordinate system (vectorized version)
v1 = B - A;
v1 = normr(v1);
v2 = cross(v1, C - A, 2);
v2 = normr(v2);
v3 = cross(v2, v1, 2);

%Projection at each time frame
for i=1:size(P,1)

    projmat = [v1(i,:)' v2(i,:)' v3(i,:)'];
    
    if(mode == 0) %normal projection
        proj.data(i,:) = (P(i,:)-A(i,:))*projmat;
    end
    if(mode == 1) %inverse projection
        if sum(sum(isnan(projmat)))==0 %no nan data
            invprojmat = pinv(projmat);
            proj.data(i,:) = P(i,:)*invprojmat + A(i,:);
        end
    end
end
proj.data(nans)=nan;
end