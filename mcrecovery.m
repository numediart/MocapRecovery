function [filledsequence] = mcrecovery(mysequence, options)

%{

Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Inputs :
- mysequence: MoCap data structure from the MoCap Toolbox
- options: optional structure, containing various optional fields:

 *binary (0 = no, any other value = yes):
    - verbose: print/display debug information. (default: 0)
    - quiet: nothing is printed/displayed (override verbose). (default: 0)
    - method1 to method5: specify any individual recovery method to use.
    (default: all to 1). As soon as a methodx parameter is specified in the
    options, other methods are set to 0.

        o method1: local interpolation
        o method2: local polynomial regression
        o method3: local GRNN
        o method4: global weighted linear regression
        o method5: weighted PCA (Gloersen et al. 2016 PLoS One)

    - spaceconstraint: use space constraint. (default: 1)
    - tripleconstraint: use space constraint with 3 neighbours. If 0, only
    the distance with the first neighbour is considered. (default: 1)
    - timeconstraint: use time constraint. (default: 1)
    - recursivefilling: each filled trajectory can be used for 
    reconstruction of the next trajectory to fill (markers with smallest 
    gaps are reconstructed first). (default: 1)
    - filtering: filter input data before process, and filter and recovered data
    (median and average symmetric sliding adaptive windows). (default: 1)
    - filterref: filter local references. (default: 1)
    - saveastsv: save recovered sequence as tsv file. (default: 0)
    - advancedordering: (under development) order possible references both
    considering maximal windowed standard deviation, and maximal distance (default:
    0)
    - advancedorderingweight: (under development) weight to ponderate
    either windowed standard deviation, or maximal distance criteria for
    ordering references (default: 20). Not used if advancedordering = 0.
    The bigger it is, the more we consider maximal windowed standard deviation over
    maximal distance

 *other:
    - fps: use this parameter if no fps is indicated in the MoCap data
    structure (mysequence). (default: 30)
    - window: filters sliding window size. Note that
    full-body motion useful information is generally below 10-20Hz. Finer
    movement may need a higher rate (e.g. finger movements). (default: 0.07)
    - presenceMin [0 - 100]: markers less present than presenceMin (in percent of 
    available frames) are removed from the sequence before the process.
    (default: 30)
    - threhsold: threshold on the distance variation (standard deviation) between the
    marker to fill and a potential reference. If the standard deviation >
    threshold, the reference is not valid for reconstruction. It can be
    used to avoid unreliable recovery. (default: inf)
    - combined_threshold: threshold on the summed distance variation (standard deviation) 
    between the marker to fill and three potential references. If the sum of the standard
    deviations > combined_threshold, the references are not valid together for
    reconstruction. (default: inf)
    - smooth: smoothing parameter for grnn regression (default: 0.3)

Output: recovered sequence

%}

%%% Parameters %%%
if nargin < 2
    options = [];
end

%if verbose = 1, debug information is printed/displayed on the process.
% (default: 0)
if isfield(options,'verbose') && isnumeric(options.verbose)
    verbose = options.verbose;
else
    verbose = 0;
end

%if quiet = 1, nothing is printed
% (default: 0)
if isfield(options,'quiet') && isnumeric(options.quiet)
    quiet = options.quiet;
else
    quiet = 0;
end

if quiet
    verbose = 0;
end

% methods to use for gap recovery (put 0 to avoid the method). 
%By default, we use all methods
% method1 : local interpolation
if isfield(options,'method1') && isnumeric(options.method1)
    method1 = options.method1;
else
    method1 = 0;
end
% method2 : local polynomial regression
if isfield(options,'method2') && isnumeric(options.method2)
    method2 = options.method2;
else
    method2 = 0;
end
% method3 : local GRNN
if isfield(options,'method3') && isnumeric(options.method3)
    method3 = options.method3;
else
    method3 = 0;
end
% method4 : global weighted linear regression
if isfield(options,'method4') && isnumeric(options.method4)
    method4 = options.method4;
else
    method4 = 0;
end
% method5 : Gloersen (PCA-based method)
if isfield(options,'method5') && isnumeric(options.method5)
    method5 = options.method5;
else
    method5 = 0;
end

% spaceconstraint : if ~= 0, soft constraint on distance of marker to fill with
% closest neighbours
% (default: 1)
if isfield(options,'spaceconstraint') && isnumeric(options.spaceconstraint)
    spaceconstraint = options.spaceconstraint;
else
    spaceconstraint = 1;
end

% tripleconstraint : if ~= 0, soft constraint on distance of marker to fill with
% 3 closest neighbours (otherwise 1 neighbour only). 
% (default: 1)
if isfield(options,'tripleconstraint') && isnumeric(options.tripleconstraint)
    tripleconstraint = options.tripleconstraint;
else
    tripleconstraint = 1;
end

% timeconstraint : if ~= 0, linear ramp correction to force trajectory
% time continuity
% (default: 1)
if isfield(options,'timeconstraint') && isnumeric(options.timeconstraint)
    timeconstraint = options.timeconstraint;
else
    timeconstraint = 1;
end

%recursivefilling : if ~= 0, each filled trajectory can be used for
%reconstruction of the next trajectory to fill (markers with smallest gaps
%are reconstructed first)
% (default: 1)
if isfield(options,'recursivefilling') && isnumeric(options.recursivefilling)
    recursivefilling = options.recursivefilling;
else
    recursivefilling = 1;
end

if isfield(options,'advancedordering') && isnumeric(options.advancedordering)
    advancedordering = options.advancedordering;
else
    advancedordering = 0;
end

if isfield(options,'advancedorderingweight') && isnumeric(options.advancedorderingweight)
    advancedorderingweight = options.advancedorderingweight;
else
    advancedorderingweight = 20;
end

%if no method is selected, use all of them (default)
if method1+method2+method3+method4+method5 == 0
    if verbose
        warning('No method indicated in sequence or in options. Default value is taken instead.');
    end
    method1=1;
    method2=1;
    method3=1;
    method4=1;
    method5=1;
end


%filter the sequence before filling (median and average symmetric sliding windows)
if isfield(options,'filtering')
    filtering = options.filtering;
else
    if verbose
        warning('No filtering parameter indicated in sequence or in options. Default value is taken instead.');
    end
    filtering = 1;
end

%filter references or not for local reconstruction (median and average symmetric sliding windows)
if isfield(options,'filterref')
    filterref = options.filterref;
else
    if verbose
        warning('No filterref parameter indicated in sequence or in options. Default value is taken instead.');
    end
    filterref = 1;
end

if isfield(mysequence,'freq')
    fps = mysequence.freq;
elseif isfield(options,'fps') && isnumeric(options.fps)
    fps = options.fps;
else
    if verbose && (filtering || filterref)
        warning('No fps indicated in sequence or in options (it is needed for filtering). Default value is taken instead.');
    end
    fps = 30;
end


if isfield(options,'window')
    window = options.window;
else
    if verbose && (filtering || filterref)
        warning('No window parameter indicated in sequence or in options (it is needed for filtering). Default value is taken instead.');
    end
    window = 0.07;
end

% markers less present than presenceMin (in % of present frames) are removed
% from the sequence before the process
if isfield(options,'presenceMin') && isnumeric(options.presenceMin)
    presenceMin = options.presenceMin;
else
    presenceMin = 30;
end
if presenceMin < 0
    warning('Invalid presenceMin parameter');
    presenceMin = 0;
elseif presenceMin >= 100
    warning('Invalid presenceMin parameter');
   presenceMin = 0; 
end

% threshold on the distance variation (standard deviation) between the
% marker to fill and a potential reference. If the standard deviation >
% threshold, the reference is not valid for reconstruction. (default: inf)
if isfield(options,'threshold') && isnumeric(options.threshold)
    threshold = options.threshold;
else
    threshold = inf;
end

% threshold on the summed distance variation (standard deviation) between the
% marker to fill and three potential references. If the sum of the standard
% deviations > combined_threshold, the references are not valid together for
% reconstruction. (default: inf)
if isfield(options,'combined_threshold') && isnumeric(options.combined_threshold)
    combined_threshold = options.combined_threshold;
else
    combined_threshold = inf;
end

%if saveastsv = 1, the recovered motion is saved in a .tsv file.
% (default: 0)
if isfield(options,'saveastsv') && isnumeric(options.saveastsv)
    saveastsv = options.saveastsv;
else
    saveastsv = 0;
end

%smoothing parameter for grnn method (default: 0.3)
if isfield(options,'smooth') && isnumeric(options.smooth)
    smooth = options.smooth;
else
    smooth = 0.3;
end

%% First check if there are missing data
[~, ~, mgrid] = mcmissing(mysequence);
if mean(mean(mgrid))==0
    %No missing data!
    filledsequence = mysequence;
    if ~quiet
        disp('There are no missing data in this sequence');
    end
    if saveastsv
        mcwritetsv2(filledsequence);
    end
    return;
end

%% Remove markers with too many missing frames

absence=mean(mgrid);
presence=(1-absence)*100;

idx = 1:mysequence.nMarkers;
markers = idx(presence>=presenceMin);
if ~quiet
    fprintf('Minimum presence: %.1f%%\n', min(presence));
    if size(markers) < mysequence.nMarkers
        fprintf('----WARNING----: Some markers have less than %.1f%% presence and are thus removed from the sequence!\n', presenceMin);
    end
end
initialsequence = mysequence;
mysequence = mcgetmarker(mysequence,markers);

%% Initialize process
mysequence.nFrames=size(mysequence.data,1);
missing = mcviewmissing(mysequence,verbose);%Diagnosis
presence=presence(markers);
name = mysequence.markerName;
x = mysequence.data;
nframes = size(x,1);
nmarkers = mysequence.nMarkers;

%% Filter data
unfiltered = mysequence;
if filtering
%     for m = 1:mysequence.nMarkers
%         mysequence.data(:,(3*m-2):3*m) = symfilter(mysequence.data(:,(3*m-2):3*m),ceil(window*fps));
%     end
    mysequence.data = symfilter(mysequence.data,ceil(window*fps));
end

%% Subtract average marker trajectory (to obtain a coordinate system moving with the subject)

cols_with_nans = any(isnan(x),1);
x_nonnan = x;
x_nonnan(:,cols_with_nans)=[];
%Center of mass of markers always presents
com.x = mean(x_nonnan(:,1:3:end),2);
com.y = mean(x_nonnan(:,2:3:end),2);
com.z = mean(x_nonnan(:,3:3:end),2);
x(:,1:3:end) = bsxfun(@minus,x(:,1:3:end),com.x);%x(:,1:3:end) - com.x %(for newer versions)
x(:,2:3:end) = bsxfun(@minus,x(:,2:3:end),com.y);
x(:,3:3:end) = bsxfun(@minus,x(:,3:3:end),com.z);
mysequence.data = x;


%% Compute 3D distances variations (standard deviation)


%Calculate distance on downsampled data for speed considerations (let's take maximum
%5000 frames randomly in data)
subsampled = mysequence;
%Keep only frames without nans
tmpdata = mysequence.data;
nanframes = any(isnan(tmpdata),2);
if size(tmpdata,1) - sum(nanframes) > 300 %Verify if we have still a sufficient number of frames
    tmpdata(nanframes,:) = [];
end
NF = min(size(tmpdata,1),5000);
randframes = randperm(size(tmpdata,1));
randframes=randframes(1:NF);
subsampled.data=tmpdata(randframes,:);
subsampled.nFrames=size(subsampled.data,1);
distances = inf(nmarkers,nmarkers);%inf avoids to reconstruct a marker with itself as reference
for i=1:nmarkers
    for j=(i+1):nmarkers
        distance = mcmarkerdist(subsampled,i,j);
        if ~advancedordering
            distances(i,j) = nanstd(distance);        
        else
            %maximal distance
            maxdist = max(distance);
            %maximal windowed standard deviation (to consider the std only
            %during movement)
            stddist = max(movstd(distance,min(250,NF/2),'omitnan'));
            %manual weighted combination of both criteria...            
            distances(i,j) = maxdist + advancedorderingweight*stddist;
        end
        distances(j,i) = distances(i,j);
    end
end

if verbose
    figure;
    imagesc(distances);colorbar();title('Standard deviation of distance between markers');
end

[dev,refs] = sort(distances);
if verbose
    figure;
    imagesc(dev);colorbar();title('Standard deviation of distance between markers (sorted)');
end

% refs are the markers that are the most linked to the marker to
% reconstruct (based on sorted distances stds for each marker).
% They are linked if the standard deviation of their distance if low (less variable
% distance). Two markers on the same segment will have a low distance std,
% but a marker on a foot and a marker on a hand will generally have a very large
% distance std.

%% Fill gaps

%Only process markers with gaps (nans)
miss = mcnumbers(name,missing)';

%Order then by presence before gap-filling (fill most present first)
[~,order] = sort(presence(miss),'descend');
miss=miss(order);

xfilled = x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gloersen (PCA-based method)
if method5
    try
        gloersen = PredictMissingMarkers_Edited(x,'Algorithm',2);
    catch
        warning('PCA-Gloersen method (5) failed and will be ignored for reconstruction');
        gloersen = x;
        if method1+method2+method3+method4 == 0
            warning('No other method than PCA-Gloersen (5) was selected. Reconstruction is thus not possible');
            warning('Consider using another method...');
            filledsequence = mysequence;
            return;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fill each missing marker
for k=miss
    
    %Extract signal
    P0=mcgetmarker(mysequence,k);
    P0filled = P0;
    linregfilled = P0;
    gloersenfilled = P0;
    p0 = P0.data;
    p0filled = P0filled.data;
    
    %Indexes of valid references for reconstruction
    validref=refs(1:end-1,k); %last one is the signal to reconstruct itself
    %Corresponding distance deviations of the valid references
    deviations = dev(1:end-1,k);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gloersen - continuity correction
    if method5
        gloersenfilled.data = gloersen(:,3*k-2:3*k);
        
        %continuity correction
        absence = isnan(P0.data(:,1))';
        %Search gaps
        begs=strfind(absence,[0 1])+1; %gaps first frames
        if absence(1)==1
            begs = [1 begs]; %in case gap at the beginning of the sequence
        end
        lasts=strfind(absence,[1 0]); %gaps last frames
        if absence(end)==1
            lasts = [lasts size(absence,2)]; %in case gap at the end on the sequence
        end
        
        lborders=begs-1; %gaps left borders
        lborders(lborders==0)=1; %in case gap at the beginning of the sequence
        rborders=lasts+1; %gaps right borders
        rborders(rborders>nframes)=P0.nFrames; %in case gap at the end on the sequence
        
        ngaps = size(lasts,2);
        temp = gloersenfilled.data;
        %Continuity correction for each gap
        if timeconstraint
            for g=1:ngaps
                id0=lborders(g); %left border
                id1=rborders(g); %right border
                gapsize = id1-id0-1;
                t1=lborders(g)+1;
                t2=rborders(g)-1;
                for ax=1:3
                    if id0 == 1 %left extremity
                        Dt1 = 0; %no correction
                    else
                        Dt1 = temp(t1,ax)-p0(t1-1,ax);
                    end
                    if id1 == nframes
                        Dt2 = 0; %no correction
                    else
                        Dt2 = temp(t2,ax)-p0(t2+1,ax);
                    end
                    if id0 == 1 && id1 ~= nframes
                        Dt1 = Dt2;
                    elseif id1 == nframes && id0 ~= 1
                        Dt2 = Dt1;
                    end
                    correction = -linspace(Dt1,Dt2,gapsize+2)'; %linear ramp correction
                    temp(t1:t2,ax) = temp(t1:t2,ax)+correction(2:end-1);
                end
            end
            gloersenfilled.data = temp;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Global linear regression (on all present markers)
    if method4
        while true %break when trajectory is filled
            
            absence = isnan(linregfilled.data(:,1))';
            if mean(absence)==1 %empty marker trajectory, no filling possible
                break;
            end
            if sum(absence)==0 %trajectory is filled!
                break;
            end
            
            %Search gaps to fill with interpolated values
            begs=strfind(absence,[0 1])+1; %gaps first frames
            if absence(1)==1
                begs = [1 begs]; %in case gap at the beginning of the sequence
            end
            lasts=strfind(absence,[1 0]); %gaps last frames
            if absence(end)==1
                lasts = [lasts size(absence,2)]; %in case gap at the end on the sequence
            end
            
            lborders=begs-1; %gaps left borders
            lborders(lborders==0)=1; %in case gap at the beginning of the sequence
            rborders=lasts+1; %gaps right borders
            rborders(rborders>nframes)=P0.nFrames; %in case gap at the end on the sequence
            
            %For each gap, check number of other markers that are present
            ngaps = size(lasts,2);
            presentDuringGaps=cell(ngaps,1);
            npresent = zeros(ngaps,1);
            for j=1:ngaps
                
                id0=lborders(j); %left border
                id1=rborders(j); %right border
                gapsequence = mctrim(mysequence,id0,id1,'frame');
                tmp = gapsequence.data(:,1:3:end);
                tmp = ~isnan(tmp);
                presentDuringGaps{j} = mean(tmp)==1;
                npresent(j) = sum(presentDuringGaps{j});
            end
            
            %Fill gaps with the most other present markers to recover them
            [nmarkers,gapidx]=max(npresent);
            if nmarkers <= 3 %not enough markers for regression
                break;
            end
            %Check if these markers are valid to fill several gaps
            refset = presentDuringGaps{gapidx};
            gapstofill = [gapidx];
            for j=1:ngaps
                if j ~= gapidx && isequal(refset, presentDuringGaps{j})
                    gapstofill = [gapstofill j];
                end
            end
            
            %Regression
            refidx = find(refset);
            tmp = mcgetmarker(mysequence,refidx);
            X = tmp.data;
            linpredictor = [ ones(nframes,1) featureNormalize(X)];
            
            %Ponderate predictors according to distance variability
            weights = distances(k,refset);
            weightScale = 20;% [mm], the sigma of the gaussian curve
            weights = exp(-weights.^2./(2*weightScale.^2));
            if max(weights) < 10^-1
                weights = ones(size(weights));
            end
            weights = [1 reshape([1 1 1]'*weights,1,[])];
            avoidids = weights <= 0.05;
            linpredictor(:,avoidids)=[];
            weights(avoidids) = [];
            linpredictor = bsxfun(@times,linpredictor,weights);
            
            regressand = p0;
            
            %remove nan frames for training
            nonnans = ~(isnan(regressand(:,1)) | any(isnan(linpredictor),2));
            regressand = regressand(nonnans,:);
            nonnanlinpredictor = linpredictor(nonnans,:);
            nonnanframes = size(regressand,1);
            
            % nmaxframes frames should be sufficient to train marker movement model
            nmaxframes = 500;
            if nonnanframes < nmaxframes
                ids=1:nonnanframes;
            elseif nonnanframes >= nmaxframes && nonnanframes < 2*nmaxframes
                %Take random frames sample to accelerate process
                ids = randperm(nonnanframes);
                ids=ids(1:nmaxframes);
                ids=sort(ids);
            elseif nonnanframes >= 2*nmaxframes
                %Downsample to accelerate process
                down = floor(nonnanframes/nmaxframes);
                ids=1:down:nonnanframes;
            end
            
            %%% Sampled data for training (faster)
            %Predictors
            predictors = nonnanlinpredictor(ids,:);%subsample
            %Regressands
            regressands = regressand(ids,:);
            [regressands, mu, sigma] = featureNormalize(regressands);
            
            p0pred=zeros(size(p0));
            for ax=1:3
                regressand = regressands(:,ax);
                warning('off');
                coeffs = regress(regressand,predictors);
                warning('on');
                pred = linpredictor*coeffs;
                pred=pred*sigma(ax)+mu(ax);
                p0pred(:,ax)=pred;
                %Continuity correction
                if timeconstraint
                    for g=gapstofill
                        t1=lborders(g)+1;
                        t2=rborders(g)-1;
                        id0=lborders(g);
                        id1 = rborders(g);
                        gapsize = id1-id0-1;
                        if id0 == 1 %left extremity
                            Dt1 = 0; %no correction
                        else
                            Dt1 = p0pred(t1-1,ax)-p0(t1-1,ax);
                        end
                        if id1 == nframes
                            Dt2=0; %no correction
                        else
                            Dt2 = p0pred(t2+1,ax)-p0(t2+1,ax);
                        end
                        if id0 == 1 && id1 ~= nframes
                            Dt1 = Dt2;
                        elseif id1 == nframes && id0 ~= 1
                            Dt2 = Dt1;
                        end
                        correction = -linspace(Dt1,Dt2,gapsize+2)'; %linear ramp correction
                        p0pred(id0:id1,ax) = p0pred(id0:id1,ax)+correction;
                    end
                end
            end
            for g=gapstofill
                %fill gap
                id0 = lborders(g);
                id1=rborders(g);
                linregfilled.data(id0:id1,:) = p0pred(id0:id1,:);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Local projection based methods
    
    %Fill with each reference UNTIL completely filled (see break in the loop)
    while(size(validref,1) >= 3) %We need at least 3 points for local projection
        
        if any(isnan(p0filled(:)))==0
            
            %Completely filled!
            xfilled(:,3*k-2:3*k) = p0filled;
            break;
        end
        
        P0filled.data = p0filled;
        %Search gaps to fill with interpolated values
        absence = isnan(p0filled(:,1))';
        if mean(absence)==1 %empty marker trajectory, no filling possible
            validref=[];
        end
        begs=strfind(absence,[0 1])+1; %gaps first frames
        if absence(1)==1
            begs = [1 begs]; %in case gap at the beginning of the sequence
        end
        lasts=strfind(absence,[1 0]); %gaps last frames
        if absence(end)==1
            lasts = [lasts size(absence,2)]; %in case gap at the end on the sequence
        end
        
        lborders=begs-1; %gaps left borders
        lborders(lborders==0)=1; %in case gap at the beginning of the sequence
        rborders=lasts+1; %gaps right borders
        rborders(rborders>nframes)=P0.nFrames; %in case gap at the end on the sequence
        
        %Look for valid references for local methods
        %References are valid if they are present during and around the gaps
        %of the signal to reconstruct, and if their distance with the signal is
        %almost constant (distance variability < threshold)
        
        %Check reference individual validity
        j=0;
        for loop=1:size(validref,1)
            
            j=j+1;
            if j > size(validref,1)
                break;
            end
            
            %Check reference presence during and around signal gap
            %If they are not present around (on borders), interpolation is
            %not possible, and if they are not present at all during, projection
            %is not possible (all nans again and we have an infinite loop)
            reftmp=mcgetmarker(mysequence,validref(j));
            nans=isnan(reftmp.data(:,1))';
            
            absentAroundGaps = nans(lborders) | nans(rborders); %at least one missing data at gap borders (no interpolation possible)
            absentDuringGaps=zeros(size(absentAroundGaps));
            for i=1:length(begs)
                absentDuringGaps(i) = mean(nans(begs(i):lasts(i)))>=1; %no data during gap
            end
            if mean(absentAroundGaps | absentDuringGaps) >= 1 || deviations(j) > threshold
                if verbose
                    fprintf('Reference marker %d not valid for reconstruction of marker %d.\n', validref(j), k);
                end
                validref(j)=[]; %remove this marker from valid references
                deviations(j)=[];
                j=j-1;
            end
        end
        
        if size(validref,1) < 3
            if ~quiet
                fprintf('Not enough valid references. Marker %d could not be completely filled.\n', k);
            end
            xfilled(:,3*k-2:3*k) = p0filled;
            break;
        end
        
        
        %Check combined validity (the 3 references must be present at the
        %same time, during and around gap, and the global deviation must be small)
        
        refselected = false;
        breaking = false;
        
        %Counters storing combinations of references
        %We must re-initialize the counters as validref may have changed.
        %It should work as previous valid combinations are no
        %more valid once used for interpolation
        c1=1;
        c2=2;
        c3=3;
        
        while refselected == false
            
            %Combined presence (0 if present, 1 if absent)
            reftmp=mcgetmarker(mysequence,validref(c1));
            nans1=isnan(reftmp.data(:,1))';
            reftmp=mcgetmarker(mysequence,validref(c2));
            nans2=isnan(reftmp.data(:,1))';
            reftmp=mcgetmarker(mysequence,validref(c3));
            nans3=isnan(reftmp.data(:,1))';
            totnans = nans1 | nans2 | nans3;
            
            %Check presence during and around gaps for the combination of
            %the three references. They must be valid together at least for
            %one gap to be selected
            absentAroundGaps = totnans(lborders) | totnans(rborders); %at least one missing data at gap borders (no interpolation possible)
            absentDuringGaps=zeros(size(absentAroundGaps));
            for i=1:length(begs)
                absentDuringGaps(i) = mean(totnans(begs(i):lasts(i)))>=1; %no data during gap
            end
            
            totdeviation = deviations(c1)+deviations(c2)+deviations(c3);
            totvaliditypergap = absentAroundGaps | absentDuringGaps;
            
            %We need at least one gap valid (at least one 0 in totvaliditypergap) for the three references at the
            %same time
            if mean(totvaliditypergap) >= 1 || totdeviation > combined_threshold
                
                %This combination of references is not valid
                if verbose
                    fprintf('Not valid references : %d, %d, %d\n', c1,c2,c3);
                end
                
                %Try next combination
                c3=c3+1;
                
                if c3 > size(validref,1) || deviations(c1)+deviations(c2)+deviations(c3) > combined_threshold
                    c2=c2+1;
                    c3=c2+1;
                end
                if c3 > size(validref,1) || deviations(c1)+deviations(c2)+deviations(c3) > combined_threshold
                    c1=c1+1;
                    c2=c1+1;
                    c3=c2+1;
                end
                if c3 > size(validref,1) || deviations(c1)+deviations(c2)+deviations(c3) > combined_threshold
                    % All possible combination tried, exit loop
                    if ~quiet
                        fprintf('Not enough valid references. Marker %d could not be completely filled.\n', k);
                    end
                    xfilled(:,3*k-2:3*k) = p0filled;
                    breaking = true;
                    refselected = true;
                end
                
            else
                %This combination of references is valid, exit loop
                refselected = true;
            end
        end
        
        if breaking == true
            % All possible combination tried
            break;
        end
        
        %Extract references
        r1=validref(c1);
        r2=validref(c2);
        r3=validref(c3);
        P1 = mcgetmarker(mysequence,r1);
        P2 = mcgetmarker(mysequence,r2);
        P3 = mcgetmarker(mysequence,r3);
        p1 = P1.data;
        p2 = P2.data;
        p3 = P3.data;
        
        %Filter data to avoid noise propagation
        if filterref
            p0 = symfilter(p0,ceil(window*fps));
            p1 = symfilter(p1,ceil(window*fps));
            p2 = symfilter(p2,ceil(window*fps));
            p3 = symfilter(p3,ceil(window*fps));
        end
        
        d1 = sqrt(sum((p1-p0).^2,2));
        d2 = sqrt(sum((p2-p0).^2,2));
        d3 = sqrt(sum((p3-p0).^2,2));
        
        if ~quiet
        fprintf('Reconstructing trajectory %s with references %s, %s and %s\n', name{k}, name{r1}, name{r2}, name{r3});
        end
        
        p0p = projection3D(p1, p2, p3, p0filled);
        
        %%% Train regression models
        
        %%% Predictors: local coordinates of P1,P2,P3, and distances
        d12 = sqrt(sum((p1-p2).^2,2));
        d13 = sqrt(sum((p1-p3).^2,2));
        d23 = sqrt(sum((p2-p3).^2,2));
        %projection of P1 is (zeros, zeros, zeros) (new origin)
        %projection of P2 is (d12, zeros, zeros)
        %projection of P3 is (P3p_x, 0, P3p_z), where :
        %d13² = P3p_x²+P3p_z²
        %d23² = P3p_z²+(P3p_x-d12)²
        p3p = projection3D(p1, p2, p3, p3);
        X=[d12 p3p(:,[1 3]) d13 d23];%
        linpredictor = [ ones(nframes,1) featureNormalize(X)];
        X=mapFeatureMat(X,2); % polynomial version
        polypredictor = [ ones(nframes,1) featureNormalize(X)];
        
        regressand = p0p;
        
        %remove nan frames for training
        nonnans = ~(isnan(regressand(:,1)) | any(isnan(linpredictor),2));
        regressand = regressand(nonnans,:);
        nonnanlinpredictor = linpredictor(nonnans,:);
        nonnanpolypredictor = polypredictor(nonnans,:);
        nonnanframes = size(regressand,1);
        
        % 2000 frames should be sufficient to train marker movement model
        nmaxframes = 1000;
        if nonnanframes < nmaxframes
            ids=1:nonnanframes;
        elseif nonnanframes >= nmaxframes && nonnanframes < 2*nmaxframes
            %Take random frames sample to accelerate process
            ids = randperm(nonnanframes);
            ids=ids(1:nmaxframes);
            ids=sort(ids);
        elseif nonnanframes >= 2*nmaxframes
            %Downsample regressand to accelerate process
            down = floor(nonnanframes/nmaxframes);
            ids=1:down:nonnanframes;
        end
        
        %%% Sampled data for training (faster)
        %Predictors
        linpredictors = nonnanlinpredictor(ids,:);%subsample
        polypredictors = nonnanpolypredictor(ids,:);%subsample
        %Regressands
        regressands = regressand(ids,:);
        [regressands, mu, sigma] = featureNormalize(regressands);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Polynomial regression on local data (3 references) with
        % continuity correction - Training
        % Used input data: local coordinates of references, and
        % distances between them (i.e. 5 non-zero components)
        coeffs = cell(3,1);
        if method2
            for ax=1:3 %regression on each axis
                regressand = regressands(:,ax);
                warning('off');
                coeffs{ax} = regress(regressand,polypredictors);
                warning('on');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NN regression on local data (3 references) with
        % continuity correction - Training
        if method3
            net = newgrnn(linpredictors',regressands',smooth);
        end
        
        
        %%% Estimate probability density function of d1, d2 and d3
        %(200 points should be enough)
        nmaxframes = 500;
        nonnans = ~isnan(d1);
        d1bis = d1(nonnans);
        if length(d1bis) > nmaxframes
            d1bis = d1bis(randperm(length(d1bis),nmaxframes));
        end
        nonnans = ~isnan(d2);
        d2bis = d2(nonnans);
        if length(d2bis) > nmaxframes
            d2bis = d2bis(randperm(length(d2bis),nmaxframes));
        end
        nonnans = ~isnan(d3);
        d3bis = d3(nonnans);
        if length(d3bis) > nmaxframes
            d3bis = d3bis(randperm(length(d3bis),nmaxframes));
        end
        
%         [pdf1,pdfi1] = ksdensity(d1bis,'NumPoints',300);
%         [pdf2,pdfi2] = ksdensity(d2bis,'NumPoints',300);
%         [pdf3,pdfi3] = ksdensity(d3bis,'NumPoints',300);
        
        
        %Fill each gap
        for j=1:size(lasts,2)
            
            if totvaliditypergap(j)==1 %(valid if 0 here)
                %We cannot fill this gap with these references
            else
                
                id0=lborders(j); %left border
                id1=rborders(j); %right border
                gapsize = id1-id0-1;
                
                %Extract refs around gaps (just before to just after)
                P1m = mctrim(P1, id0, id1,'frame');
                P2m = mctrim(P2, id0, id1,'frame');
                P3m = mctrim(P3, id0, id1,'frame');
                P0m = mctrim(P0filled ,id0, id1,'frame');
                p1m = P1m.data;
                p2m = P2m.data;
                p3m = P3m.data;
                p0m = P0m.data;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Local interpolation on 3 references
                if method1
                    invproj = local_interpolation(P1m,P2m,P3m,P0m);
                    locinterp=invproj.data;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Local polynomial regression - prediction and continuity correction
                if method2
                    temp = p0p(id0:id1,:);
                    for ax=1:3
                        pred = polypredictor(id0:id1,:)*coeffs{ax};
                        pred=pred*sigma(ax)+mu(ax);
                        temp(:,ax)=pred;
                        if timeconstraint
                            if id0==1
                                Dt1=0;
                            else
                                Dt1 = temp(1,ax)-p0p(id0,ax);
                            end
                            if id1 == nframes
                                Dt2=0;
                            else
                                Dt2 = temp(end,ax)-p0p(id1,ax);
                            end
                            
                            if id0 == 1 && id1 ~= nframes
                                Dt1 = Dt2;
                            elseif id1 == nframes && id0 ~= 1
                                Dt2 = Dt1;
                            end
                            
                            correction = -linspace(Dt1,Dt2,gapsize+2)'; %linear ramp correction
                            temp(:,ax) = temp(:,ax)+correction;
                        end
                    end
                    %Inverse projection
                    locreg = projection3D(p1m, p2m, p3m, temp, 1);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Local NN regression - prediction and continuity correction
                if method3
                    %temp = p0plocNN(id0:id1,:);
                    tmppredictor = linpredictor(id0:id1,:);
                    pred = sim(net,tmppredictor');
                    pred=pred';
                    pred=bsxfun(@times, pred, sigma);
                    temp=bsxfun(@plus, pred, mu);
                    if timeconstraint
                        for ax=1:3
                            if id0==1
                                Dt1=0;
                            else
                                Dt1 = temp(1,ax)-p0p(id0,ax);
                            end
                            if id1 == nframes
                                Dt2=0;
                            else
                                Dt2 = temp(end,ax)-p0p(id1,ax);
                            end
                            
                            if id0 == 1 && id1 ~= nframes
                                Dt1 = Dt2;
                            elseif id1 == nframes && id0 ~= 1
                                Dt2 = Dt1;
                            end
                            
                            correction = -linspace(Dt1,Dt2,gapsize+2)'; %linear ramp correction
                            temp(:,ax) = temp(:,ax)+correction;
                        end
                    end
                    %Inverse projection
                    locNN = projection3D(p1m, p2m, p3m, temp, 1);
                    
%                     %%% temporary code for time constraint effect visualization
%                     origdata = p0p;
%                     %origdata = origdata(:,3);
%                     corrected1 = origdata;
%                     notimeconst = bsxfun(@plus, pred, mu);
%                     corrected1(id0+1:id1-1,:) = notimeconst(2:end-1,:);
%                     corrected2 = origdata;
%                     corrected2(id0+1:id1-1,:) = temp(2:end-1,:);
%                     figure;
%                     tmprange = 1:length(origdata);
%                     ax = 3;
%                     plot(corrected1(tmprange,ax));hold on;plot(corrected2(tmprange,ax));plot(origdata(tmprange,ax));
%                     figure;
%                     for ax = 1:3
%                         subplot(3,1,ax);
%                         plot(corrected1(tmprange,ax));hold on;plot(corrected2(tmprange,ax));plot(origdata(tmprange,ax));
%                     end
                    
                end
                
                
                
%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Kernel smoothing regression - prediction and continuity correction
%                 if method6
%                     temp = p0p(id0:id1,:);
%                     for ax=1:3
%                         %pred = polypredictor(id0:id1,:)*coeffs{ax};
%                         %pred = kdtree_ksrmv(linpredictors,regressands(:,ax),[],linpredictor(id0:id1,:));
%                         pred = ksrmv(linpredictors,regressands(:,ax),[],linpredictor(id0:id1,:));
%                         pred=pred.f*sigma(ax)+mu(ax);
%                         temp(:,ax)=pred;
%                         if timeconstraint
%                             if id0==1
%                                 Dt1=0;
%                             else
%                                 Dt1 = temp(1,ax)-p0p(id0,ax);
%                             end
%                             if id1 == nframes
%                                 Dt2=0;
%                             else
%                                 Dt2 = temp(end,ax)-p0p(id1,ax);
%                             end
%                             
%                             if id0 == 1 && id1 ~= nframes
%                                 Dt1 = Dt2;
%                             elseif id1 == nframes && id0 ~= 1
%                                 Dt2 = Dt1;
%                             end
%                             
%                             correction = -linspace(Dt1,Dt2,gapsize+2)'; %linear ramp correction
%                             temp(:,ax) = temp(:,ax)+correction;
%                         end
%                     end
%                     %Inverse projection
%                     locksr = projection3D(p1m, p2m, p3m, temp, 1);
%                 end
%                 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Model averaging
                %(based on distance probabilities distribution estimation)
                
                %methods
                method=[];
                met=1;
                if method1 && any(any(~isnan(locinterp)))
                    method{met}=locinterp;
                    met=met+1;
                end
                if method2 && any(any(~isnan(locreg)))
                    method{met}=locreg;
                    met=met+1;
                end
                if method3 && any(any(~isnan(locNN)))
                    method{met}=locNN;
                    met=met+1;
                end
                if method4 && any(any(~isnan(linregfilled.data(begs(j):lasts(j),:))))
                    method{met}=linregfilled.data(id0:id1,:);
                    met=met+1;
                end
                if method5 && any(any(~isnan(gloersenfilled.data(begs(j):lasts(j),:))))
                    method{met}=gloersenfilled.data(id0:id1,:);
                    met=met+1;
                end
%                 if method6 && any(any(~isnan(locksr)))
%                     method{met}=locksr;
%                     met=met+1;
%                 end
                
%                 %add pchip is gap is small enough
%                 if lasts(j) - begs(j) < 0.1*fps
%                     pchipfilled = naninterp(p0,'pchip');
%                     method{met}=pchipfilled(id0:id1,:);
%                     met=met+1;
%                 end
                
                nmethods = length(method);
                weight=cell(nmethods,1);
                meanweight = zeros(nmethods,1);
                if nmethods > 1
                    for i=1:nmethods
                        
                        %comment the following line if you do not want to
                        %filter individually recovered gaps
                        method{i} = symfilter(method{i},ceil(window*fps));
                        
                        %distances of reconstructed trajectory with references during the gap
                        dist1 = sqrt(sum((method{i}-p1m).^2,2));
                        dist2 = sqrt(sum((method{i}-p2m).^2,2));
                        dist3 = sqrt(sum((method{i}-p3m).^2,2));
                        
                        %probabilities of each reconstructed trajectory (at each frame)
                        %according to its distance with each reference
%                         w1 = max(interp1(pdfi1,pdf1,dist1,'pchip',0),10^-100);%safeguard (wi>0)
%                         w2 = max(interp1(pdfi2,pdf2,dist2,'pchip',0).^(1/10),10^-100);
%                         w3 = max(interp1(pdfi3,pdf3,dist3,'pchip',0).^(1/20),10^-100);
                        w1=ksdensity(d1bis,dist1);
                        w2=ksdensity(d2bis,dist2);
                        w3=ksdensity(d3bis,dist3);
                        weight{i} = w1.*w2.*w3;
                        %w = mvksdensity([d1bis d2bis d3bis],[dist1 dist2 dist3]);
                        %weight{i} = w;
                        weight{i} = max(weight{i},10^-200);%safeguard (w>0)
                        %filter weights to avoid flickering when a bad
                        %recovery has by chance a good score
                        if verbose
                            figure(7);subplot(nmethods,1,i);plot(weight{i})
                        end
                        weight{i} = symfilter(weight{i},min(ceil(fps*window),ceil(gapsize/2)));
                        %weight{i} = movmedian(weight{i},min(ceil(fps),ceil(gapsize/2)));
                        %weight{i} = movmean(weight{i},min(ceil(fps),ceil(gapsize/2)));
                        %for a very short gap, just keep one value (avoid
                        %flickering)
                        %                         if gapsize < window*fps
                        %                            weight{i}(:) = median(weight{i});
                        %                         end
                        
                        meanweight(i) = median(weight{i});
                    end
                    
                    if verbose
                        figure(6);clf;
                        for m=1:nmethods
                            subplot(nmethods,1,m);plot(weight{m});
                        end
                    end
                    
                    temp=zeros(gapsize+2,3);
                    for g=1:(gapsize+2)
                        sumweights=0;
                        for i=1:nmethods
                            temp(g,:)=temp(g,:)+weight{i}(g)*method{i}(g,:);
                            sumweights = sumweights + weight{i}(g);
                        end
                        temp(g,:)=temp(g,:)/sumweights;
                    end
                else
                    temp = method{1};
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if spaceconstraint
                    %Add soft constraint on p1
                    %                 m1=nanmean(d1);
                    %                 s1=nanstd(d1);
                    %                 R1 = max(d1);
                    %                 r1 = min(d1);
                    %                 R1 = min(m1+2*s1,R1);
                    %                 r1 = max(m1-2*s1,r1);
                    R1=quantile(d1,0.99);
                    r1=quantile(d1,0.01);
                    %cdf1=cumsum(pdf1)/sum(pdf1);
                    %confidence interval [r1;R1]
                    %R1=pdfi1(find(cdf1<=0.99,1,'last'));
                    %r1=pdfi1(find(cdf1>=0.01,1,'first'));
                    %R1 = pdfi1(95);
                    %r1 = pdfi1(5);
                    if ~tripleconstraint
                        for g=1:(gapsize+2)
                            P0g = temp(g,:);
                            P1g=p1m(g,:);
                            %Check distance between p0 and p1; if it is outside the
                            %acceptable range ([r1,R1]), we move p0 along the line
                            %p0-p1 so that the distance comes in the range.
                            dist1g = norm(P0g-P1g);
                            if dist1g > R1
                                d=R1;
                                v1 = (P0g-P1g)/dist1g;
                                Pp = P1g+v1*d;
                            elseif dist1g < r1
                                d=r1;
                                v1 = (P0g-P1g)/dist1g;
                                Pp = P1g+v1*d;
                            else
                                Pp=P0g;
                            end
                            temp(g,:)=Pp;
                        end
                    else
                        %cdf2=cumsum(pdf2)/sum(pdf2);
                        %cdf3=cumsum(pdf3)/sum(pdf3);
                        %R2=pdfi2(find(cdf2<=0.95,1,'last'));
                        %r2=pdfi2(find(cdf2>=0.05,1,'first'));
                        %R3=pdfi3(find(cdf3<=0.95,1,'last'));
                        %r3=pdfi3(find(cdf3>=0.05,1,'first'));
                        R2=quantile(d2,0.99);
                        r2=quantile(d2,0.01);
                        R3=quantile(d3,0.99);
                        r3=quantile(d3,0.01);
                        for g=1:(gapsize+2)
                            temp(g,:) = recursivespaceconstraint(temp(g,:),p1m(g,:), r1, R1, p2m(g,:), r2, R2, p3m(g,:), r3, R3);
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Replace gap (and gap only!) with filled values
                p0filled(id0:id1,:) = temp;
                %Check distance variability for filled data
                dist1 = sqrt(sum((temp-P1m.data).^2,2));
                dist2 = sqrt(sum((temp-P2m.data).^2,2));
                deviation1 = nanstd(dist1);
                deviation2 = nanstd(dist2);
                if deviation1 > deviations(c1) || deviation2 > deviations(c2)
                    if verbose
                        disp('----WARNING----: new distance variability between marker and references is unexpectedly high. You should check the result!');
                    end
                end
            end
        end
        
        if verbose
            dime=3;
            figure(5);
            subplot(2,1,1);
            plot(p0(:,dime));title('original');
            subplot(2,1,2);
            plot(p0filled(:,dime));title('filled');
            %system('pause');
        end
        
        if(any(isnan(p0filled(:)))==0)
            
            %Completely filled!
            xfilled(:,3*k-2:3*k) = p0filled;
            
            if recursivefilling
                mysequence.data(:,3*k-2:3*k) = p0filled;
            end
            break;
        else
            if ~quiet
                disp('Not filled yet!');
            end
        end
        
    end %for each reference UNTIL completely filled
    
end %for each missing marker

%% Add com trajectory
xfilled(:,1:3:end) = bsxfun(@plus,xfilled(:,1:3:end),com.x); %xfilled(:,1:3:end) + com.x; % works for recent versions of matlab
xfilled(:,2:3:end) = bsxfun(@plus,xfilled(:,2:3:end),com.y);
xfilled(:,3:3:end) = bsxfun(@plus,xfilled(:,3:3:end),com.z);
filledsequence = mysequence;
filledsequence.data = xfilled;

%% Filter data

if filtering
%     for m = 1:filledsequence.nMarkers
%         filledsequence.data(:,(3*m-2):3*m) = symfilter(filledsequence.data(:,(3*m-2):3*m),ceil(window*fps));
%     end
    filledsequence.data = symfilter(filledsequence.data,ceil(window*fps));
end

%% Put back uncorrupted data (not filtered)
filledsequence.data(~isnan(unfiltered.data)) = unfiltered.data(~isnan(unfiltered.data));

%% Add back dismissed markers for reconstruction (presence below threshold)
%needed to the same number of markers

%comment the following line if you want to keep the dismissed markers
%(otherwise they as filled with nans for the whole sequence)
initialsequence.data = initialsequence.data*nan;


for m=1:length(markers)
    mark = markers(m);
    initialsequence.data(:,3*mark-2:3*mark) = filledsequence.data(:,3*m-2:3*m);
end

filledsequence = initialsequence;




%% Save corrected sequence

if saveastsv
    mcwritetsv2(filledsequence);
end


end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%internal functions

%recusrive projections on confidence intervals (CI) until all CI are respected
%softer constraint first
function out = recursivespaceconstraint(in, p1, r1, R1, p2, r2, R2, p3, r3, R3, maxit)

if nargin < 11
    maxit = 5;
end

out = in;

for it = 1:maxit
    
    dist1 = norm(out-p1);
    dist2 = norm(out-p2);
    dist3 = norm(out-p3);
    %end process if all CI conditions are respected
    if dist1 >= r1 && dist1 <= R1 && dist2 >= r2 && dist2 <= R2 && dist3 >= r3 && dist3 <= R3
        break;
    end
    
    %projection on p3 confidence interval
    dist3 = norm(out-p3);
    if dist3 > R3 %high bound
        d=R3;
        dir3 = (out-p3)/dist3;
        out = p3+dir3*d;
    elseif dist3 < r3 %low bound
        d=r3;
        dir3 = (out-p3)/dist3;
        out = p3+dir3*d;
    end
    %projection on p2 confidence interval
    dist2 = norm(out-p2);
    if dist2 > R2 %high bound
        d=R2;
        dir2 = (out-p2)/dist2;
        out = p2+dir2*d;
    elseif dist2 < r2 %low bound
        d=r2;
        dir2 = (out-p2)/dist2;
        out = p2+dir2*d;
    end
    %projection on p1 confidence interval
    dist1 = norm(out-p1);
    if dist1 > R1 %high bound
        d=R1;
        dir1 = (out-p1)/dist1;
        out = p1+dir1*d;
    elseif dist1 < r1 %low bound
        d=r1;
        dir1 = (out-p1)/dist1;
        out = p1+dir1*d;
    end
end

end


 
 

