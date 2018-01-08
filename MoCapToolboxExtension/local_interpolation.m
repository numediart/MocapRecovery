function invproj = local_interpolation(A, B, C, P)

%{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Interpolate the marker 3D trajectory P in a local coordinate system 
defined by A, B and C
%}


%Projection
proj = mcprojection(A, B, C, P);

%If the first frame is missing (gap at the beginning of the track),
%we put some data so we can still interpolate; the best we can do
%is to copy the first non-nan value
if isnan(P.data(1,1))
    try
        nonnan = find(isnan(proj.data(:,1))==0,1,'first');% first non-nan value
        if(isempty(nonnan))
            
        end
        proj.data(1,1) = proj.data(nonnan,1);
        proj.data(1,2) = proj.data(nonnan,2);
        proj.data(1,3) = proj.data(nonnan,3);
    catch
        disp('Error in local_interpolation (missing data at first frame)');
        disp(proj.data(1,1));
        system('pause');
    end
    
end
%Idem if the last frame is missing
if isnan(P.data(end,1))
    try
        nonnan = find(isnan(proj.data(:,1))==0,1,'last');% last non-nan value
        proj.data(end,1) = proj.data(nonnan,1);
        proj.data(end,2) = proj.data(nonnan,2);
        proj.data(end,3) = proj.data(nonnan,3);
    catch
        disp('Error in local_interpolation (missing data at last frame)');
        disp(proj.data(end,1));
        system('pause');
    end
end

%Interpolation
interp=proj;
interp.data = naninterp(proj.data,'linear');

%Inverse projection
invproj = mcprojection(A, B, C, interp, 1);

end