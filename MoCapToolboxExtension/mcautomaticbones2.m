function connections = mcautomaticbones2(mysequence,nbones,weight)

%{

Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Create automatically a connection matrix from MoCap data

%}

if nargin < 2
    nbones = round(2.5*mysequence.nMarkers); % markers distance (in mm) under which a connection is possible
end

if nargin < 3
    weight = 15; % markers distance (in mm) under which a connection is possible
end
%%% Initialize process
mysequence.nFrames=size(mysequence.data,1);
nmarkers = mysequence.nMarkers;

%%% Compute 3D distances variations (standard deviation)

connections = [];
%Calculate distance on downsampled data for speed considerations (let's take maximum
%5000 frames randomly in data)
subsampled = mysequence;
%Keep only frames without nans
data = mysequence.data;
nanframes = any(isnan(data),2);
if size(data,1) - sum(nanframes) > 300 %Verify if we have still a sufficient number of frames
    data(nanframes,:) = [];
end

NF = min(size(data,1),5000);
randframes = randperm(size(data,1));
randframes=randframes(1:NF);
subsampled.data=data(randframes,:);
subsampled.nFrames=size(subsampled.data,1);
%distances = inf(nmarkers,nmarkers);%inf avoids to reconstruct a marker with itself as reference
scores = ones(nmarkers,nmarkers)*inf;

for i=1:nmarkers
    for j=(i+1):nmarkers
        distance = mcmarkerdist(subsampled,i,j);
        %distances(i,j) = nanstd(distance);
        %distances(j,i) = distances(i,j);
        %         if nanstd(distance) < threshold
        %             connections = [connections; i j];
        %         end
        %
        %more robust: windowed standard deviation
        %         if max(movstd(distance,min(500,NF/2))) < threshold && max(distance) < maxdist
        %              connections = [connections; i j];
        %         end
        
        %weighted combination of criteria
        %maximal distance
        mymaxdist = max(distance);
        %maximal windowed standard deviation (to consider the std only
        %during movement)
        stddist = max(movstd(distance,min(250,NF/2),'omitnan'));
        %manual weighted combination of both criteria...
        scores(i,j) = mymaxdist + weight*stddist;      
    end
end

tmpscores = sort(scores(:));
threshold = tmpscores(nbones);

for i=1:nmarkers
    for j=(i+1):nmarkers
        if scores(i,j) <= threshold
            connections = [connections; i j];
        end        
    end
end

end