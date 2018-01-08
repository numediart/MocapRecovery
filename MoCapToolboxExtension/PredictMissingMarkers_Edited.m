

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Modified version of the function PredictMissingMarkers, from:
%{
 Gløersen Ø, Federolf P (2016) Predicting Missing Marker Trajectories in
 Human Motion Data Using Marker Intercorrelations. PLoS ONE 11(3):
 e0152616. https://doi.org/10.1371/journal.pone.0152616
 
 Modifications are commented in the code, and are mainly aimed at reducing the
 processing time. The original function could be used instead.
 
%}

function [GapfilledDataSet] = PredictMissingMarkers_Edited(Data_gaps,varargin)
% PredictMissingMarkers
% searches for gaps ("NaN"s) in the marker coordinates and fills these gaps
% using the correlation structure in the marker coordinates.
%
% Input:
% Required:
% Data_gaps: matrix with marker data organized in the form
%            [ x1(t1), y1(t1), z1, x2, y2, z2, x3, ..., zn(t1)
%              x1(t2), y1(t2), ...
%               ...                                     ...
%              x1(tm), y1(tm), ...                    , zn(tm)]
%
% Optional parameter - value pairs:
% 'Algorithm': Reconstruction strategy for gaps in multiple markers
%                   (Integers 1 or 2, corresponding to strategies R1 and R2)
% 'WeightScale':    Parameter "sigma" for determining weights. (Default: 200)
% 'MMweight':       Weight on missing marker (Default: 0.02)
% 'DistalThreshold':    Cutoff distance for distal markers in R2, relative to
%                   average Euclidean distance between all markers. (Default: 0.5)
% 'MinCumSV':       Cumumlative sum of normalized singular values that determines the
%                   number of PC-vectors included in the analysis (Default: 0.99)
%
%
% Output:
% GapfilledDataSet: Matrix in the same form as the input matrix
%                           where the missing data frames is replaced by
%                           reconstructed marker trajectories


%% Set program parameters:

%default values
parser = inputParser;
defaultAlgorithm = 2;
defaultWeightScale = 200; % [mm]
defaultMMweight = 0.02;
defaultDistalThreshold = 0.5;
defaultMinCumSV = 0.99;

expectedAlgorithm = [1 2];
checkRecAlg = @(x) isnumeric(x) && ismember(x,expectedAlgorithm);
checkCumVar = @(x) isnumeric(x) && x>0 && x<=1;
checkMMweight = @(x) isnumeric(x) && x>0;

addRequired(parser,'Data_gaps',@isnumeric);
addOptional(parser,'WeightScale',defaultWeightScale,@isnumeric);
addOptional(parser,'MMweight',defaultMMweight,checkMMweight);
addOptional(parser,'DistalThreshold',defaultDistalThreshold,@isnumeric);
addOptional(parser,'MinCumSV',defaultMinCumSV,checkCumVar);
addOptional(parser,'Algorithm',defaultAlgorithm,checkRecAlg);

%set parameters:
parse(parser,Data_gaps,varargin{:});

weightScale = parser.Results.WeightScale;
MMweight = parser.Results.MMweight;
DistalThreshold = parser.Results.DistalThreshold;
MinCumSV = parser.Results.MinCumSV;


%% Check input file for gaps

[~,columns] = size(Data_gaps);
% Detect which columns have gaps and where
ColumnsWithGaps = find(any(isnan(Data_gaps),1));
MarkersWithGaps = ColumnsWithGaps(3:3:end)./3;
%FramesWithGaps = find(any(isnan(Data_gaps),2));

if isempty(find(any(isnan(Data_gaps),2),1))
    warning('Submitted data does not appear to have any gaps. Make sure that gaps are represented by NaNs')
    GapfilledDataSet = Data_gaps;
    return
end

DataWithoutCompromizedMarkers = Data_gaps;
DataWithoutCompromizedMarkers(:,ColumnsWithGaps)=[];

%% Subtract average marker trajectory (to obtain a coordinate system moving with the subject)

mean_trajectory.x = ...
    mean(DataWithoutCompromizedMarkers(:,...
    1:3:columns-size(ColumnsWithGaps,2) ),2);
mean_trajectory.y = ...
    mean(DataWithoutCompromizedMarkers(:,...
    2:3:columns-size(ColumnsWithGaps,2) ),2);
mean_trajectory.z = ...
    mean(DataWithoutCompromizedMarkers(:,...
    3:3:columns-size(ColumnsWithGaps,2) ),2);

Data_gaps(:,1:3:columns) = ...
    Data_gaps(:,1:3:columns)-...
    repmat(mean_trajectory.x,1,columns/3);
Data_gaps(:,2:3:columns) = ...
    Data_gaps(:,2:3:columns)-...
    repmat(mean_trajectory.y,1,columns/3);
Data_gaps(:,3:3:columns) = ...
    Data_gaps(:,3:3:columns)-...
    repmat(mean_trajectory.z,1,columns/3);

%% Allocate space for output matrix:
GapfilledDataSet = Data_gaps;

%% Choose reconstruction strategy
switch parser.Results.Algorithm
    case 1
        % disp('R1 strategy')
        ReconstructedFullDataSet = reconstruct(Data_gaps);
        IndexesWithGaps = find(isnan(Data_gaps));
        GapfilledDataSet(IndexesWithGaps) = ReconstructedFullDataSet(IndexesWithGaps);
    case 2
        % disp('R2 strategy')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Edit from Mickaël Tits (02/2017) :
        %distance2marker can be called only once (faster)
        
        Distances = distance2marker(Data_gaps,ColumnsWithGaps); %weightvector is now the Eucl. dist from the gapped marker(s) to all other markers
        m = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        for i=MarkersWithGaps
            % remove columns distal to marker i
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Edit from Mickaël Tits (02/2017) :
            %avoid unnecessary calls of distance2marker (as it is time-consuming)
            
            %EuclDist2Markers = distance2marker(Data_gaps,3*i-2:3*i);
            m=m+1;
            EuclDist2Markers = Distances(m,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            thresh =  DistalThreshold.*mean(EuclDist2Markers);
            Cols2Zero = reshape(repmat(EuclDist2Markers,3,1),1,columns)>thresh ...
                & any(isnan(Data_gaps),1);
            Data_gaps_removedCols = Data_gaps;
            Data_gaps_removedCols(:,Cols2Zero) = 0;
            Data_gaps_removedCols(:,3*i-2:3*i) = Data_gaps(:,3*i-2:3*i);
            %Find overlapping gaps with marker i, and remove them from the matrix
            GappedFrames_i  = isnan(Data_gaps(:,3*i));
            for ii=MarkersWithGaps(MarkersWithGaps~=i)
                if any(isnan(Data_gaps_removedCols(GappedFrames_i,3*ii)))
                    Data_gaps_removedCols(:,3*ii-2:3*ii) = 0;
                end
            end % end removing overlapping gaps
            
            %           %remove frames with incomplete information:
            FramesNoGaps = find(~any(isnan(Data_gaps_removedCols),2));
            %add frames with gaps in marker i to the end of the matrix:
            frames2rec = find(any(isnan(Data_gaps(:,3*i)),2));
            CompleteAndGappedFrames = ([FramesNoGaps' frames2rec']);
            FillFrames = length(FramesNoGaps)+1:length(CompleteAndGappedFrames);
            TempReconstructedData = reconstruct(Data_gaps_removedCols(CompleteAndGappedFrames,:));
            GapfilledDataSet(frames2rec,3*i-2:3*i) = ...
                TempReconstructedData(FillFrames,3*i-2:3*i);
        end %end filling marker i
        
    otherwise
        disp('Error: invalid reconstruction approach')
        GapfilledDataSet = nan;
end

%% Add mean trajectory

GapfilledDataSet(:,1:3:columns) = ...
    GapfilledDataSet(:,1:3:columns) + ...
    repmat(mean_trajectory.x,1,columns/3);
GapfilledDataSet(:,2:3:columns) = ...
    GapfilledDataSet(:,2:3:columns) + ...
    repmat(mean_trajectory.y,1,columns/3);
GapfilledDataSet(:,3:3:columns) = ...
    GapfilledDataSet(:,3:3:columns) + ...
    repmat(mean_trajectory.z,1,columns/3);

% Done filling gaps

%% reconstruction function
    function Reconstruction = reconstruct(Data2reconstruct)
        %Determine frames and columns with gaps
        ColsWithGaps = find(any(isnan(Data2reconstruct),1));
        %find the weight vector based on Euclidean distances between markers
        weightvector = distance2marker(Data_gaps,ColsWithGaps); %weightvector is now the Eucl. dist from the gapped marker(s) to all other markers
        % for R1-approach: keep only the smallest distances for weights
        if size(weightvector,1)>1
            weightvector = min(weightvector);
        end
        weightvector = exp(-weightvector.^2./(2*weightScale.^2));
        weightvector(ColsWithGaps(3:3:end)/3) = MMweight;
        
        %define matrices needed for reconstruction:
        M = Data2reconstruct;
        
        M_zeros = M;
        M_zeros(:,ColsWithGaps) = 0;
        
        N_no_gaps = Data2reconstruct(~any(isnan(M),2),:);
        
        N_zeros = N_no_gaps;
        N_zeros(:,ColsWithGaps) = 0;
        
        %normalize to unit variance, then multiply by weighting vector:
        mean_N_no_gaps    = mean(N_no_gaps,1);
        mean_N_zeros = mean(N_zeros,1);
        stdev_N_no_gaps   = std(N_no_gaps,1,1);
        stdev_N_no_gaps(stdev_N_no_gaps==0) = 1;
        %the zeroed colums should remain zero. Set the stdev_N_no_gaps to
        % one for the zeroed columns, so that they remain zero after
        % normalization
        M_zeros = (M_zeros - repmat(mean_N_zeros,size(M_zeros,1),1))./...
            repmat(stdev_N_no_gaps,size(M_zeros,1),1).*...
            repmat(reshape([1 1 1]'*weightvector,1,[]),...
            size(M_zeros ,1),1);
        
        N_no_gaps =(N_no_gaps-repmat(mean_N_no_gaps,size(N_no_gaps ,1),1))./...
            repmat(stdev_N_no_gaps,size(N_no_gaps ,1),1).*...
            repmat(reshape([1 1 1]'*weightvector,1,[]),...
            size(N_no_gaps ,1),1);
        
        N_zeros = (N_zeros - repmat(mean_N_zeros,size(N_zeros ,1),1))./...
            repmat(stdev_N_no_gaps,size(N_no_gaps ,1),1).*...
            repmat(reshape([1 1 1]'*weightvector,1,[]),...
            size(N_no_gaps ,1),1);
        
        
        [PC_vectors_no_gaps,sqrtEigvals_no_gaps] = PCA(N_no_gaps);
        
        [PC_vectors_zeros,sqrtEigvals_zeros] = PCA(N_zeros);
        
        % Select the number of PV-vectors to include in the analysis
        n_eig = [find(cumsum(sqrtEigvals_no_gaps)>=...
            MinCumSV*sum(sqrtEigvals_no_gaps),1,'first');
            find(cumsum(sqrtEigvals_zeros)>=...
            MinCumSV*sum(sqrtEigvals_zeros),1,'first')];
        n_eig = max(n_eig);
        
        PC_vectors_no_gaps = PC_vectors_no_gaps(:,1:n_eig);
        PC_vectors_zeros = PC_vectors_zeros(:,1:n_eig);
        % Calculate Transformation Matrix
        T = PC_vectors_no_gaps'*PC_vectors_zeros;
        
        % Transform Data first into incomplete-, then into full-PC basis system.
        ReconstructedData = M_zeros*PC_vectors_zeros*T*PC_vectors_no_gaps';
        % (Equation 1 in Federolfs 2013 paper)
        
        % Reverse normalization
        ReconstructedData = repmat(mean_N_no_gaps,size(Data2reconstruct,1),1)...
            + ReconstructedData.*repmat(stdev_N_no_gaps,size(ReconstructedData ,1),1)./...
            repmat(reshape([1 1 1]'*weightvector,1,[]),size(M_zeros ,1),1);
        %prepare  output
        Reconstruction = Data2reconstruct;
        for j = ColumnsWithGaps
            Reconstruction(:,j) = ReconstructedData(:,j);
            
        end
        
    end %end function reconstruct


end %end PredictMissingMarkers

% PredictMissingMarkers subfunctions:
function [PC,sqrtEV] = PCA(Data)
[N,~] = size(Data);
Y = Data / sqrt(N-1);
[~,sqrtEV,PC] = svd(Y,'econ');
sqrtEV = diag(sqrtEV);
end

%spatial distance function
function [distArray] = distance2marker(MarkerData,colWithGaps)

[n,~] = size(MarkerData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Edit from Mickaël Tits (02/2017) : downsample for speed
nframesmax = 500;
if n < nframesmax
    ids=1:n;
elseif n >= nframesmax && n < 2*nframesmax
    %Take random frames sample to accelerate process
    ids = randperm(n);
    ids=ids(1:nframesmax);
    ids=sort(ids);
elseif n >= 2*nframesmax
    %Downsample regressand to accelerate process
    down = floor(n/nframesmax);
    ids=1:down:n;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MarkerData = MarkerData(ids,:);
[n,m] = size(MarkerData);

MarkerWithGaps = colWithGaps(3:3:end)/3;
nMarkerWithGaps = length(MarkerWithGaps);
MarkerData = reshape(MarkerData',3,m/3,n);

distArray = nan(nMarkerWithGaps,m/3,n);
for i=1:n
    distArray(:,:,i) =pdist2(MarkerData(:,MarkerWithGaps,i)',...
        MarkerData(:,:,i)','euclidean');
end
distArray = nanmean(distArray,3);

end

