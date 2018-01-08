function output = symfilter(input,windowsize,usemedian,usemean)

% This function was originally designed to filter MoCap data
% (as a matrix [nb of frames , 3*nb of markers]), while keeping nan values,
% and keeping signal edges and gap borders (near nans) at the right place.
% To that end, a sliding symmetric window is applied on each column of the
% signal, and the size of the window is shrinked near edges and gap borders
% to preserve symmetry.
% A non-symmetric window would lead to a delay in the filtered signal.
% This function filters the input data successively with a median moving
% window and a mean moving window.
% The median removes peaks (glitches and high frequency noise), while the
% mean smoothes data.

% input arguments:

% input : matrix of MoCap data [nb of frames , 3*nb of markers]

% windowsize : size in frames of the sliding window (see matlab
% documentation on movmedian and movmean for more details).
% Note: full-body human movement is generally far below 20Hz, so you can
% easily use a window of 1/20 second, i.e. a size of 1/20 * fps

% usemedian : use or not median filtering (default 1). Put it to 0 to avoid
% median filtering
% usemean : use or not mean filtering (default 1). Put it to 0 to avoid
% mean filtering

% Copyright (C) <2017)> <Mickael Tits>

% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.

% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
% OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
% CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

if nargin < 3
    usemedian = 1;
    usemean = 1;
end


transposed = 0;
if size(input,1) == 1 && size(input,2) > 1
    input = input';
    transposed = 1;
end

output = input;
if windowsize < 2
    return;
end

for col = 1:size(input,2)
    inputcol = input(:,col);
    
    if usemedian

    %%% Median filtering
    % filter of windowsize (first step)
    % option 'fill' allows to leave edges and gap borders as nan
    
        tmp = movmedian(inputcol,windowsize,'Endpoints','fill');
        
        %filter one by one the rest with an adaptative symmetric window (near edges and nans)
        rest = isnan(tmp(:,1)) & ~isnan(inputcol(:,1));
        nans = isnan(inputcol(:,1));
        restid = find(rest);
        %nans and edges
        nanids = [0; find(nans); (size(inputcol,1) + 1)];
        for i=1:length(restid)
            id=restid(i);
            %find closest nan to define largest acceptable size of symmetric window
            maxdist = min(abs(nanids-id))-1;
            if maxdist < 2
                tmp(id,:) = inputcol(id,:);
            else
                window = id-maxdist:id+maxdist;
                tmp(id,:) = median(inputcol(window,:));
            end
        end
        
        inputcol = tmp;
    end
    
    if usemean
        
        %%% Mean filtering
        % filter of windowsize (first step)
        tmp = movmean(inputcol,windowsize,'Endpoints','fill');
        
        %filter one by one the rest with an adaptative symmetric window (near edges and nans)
        rest = isnan(tmp(:,1)) & ~isnan(inputcol(:,1));
        nans = isnan(inputcol(:,1));
        restid = find(rest);
        %nans and edges
        nanids = [0; find(nans); (size(inputcol,1) + 1)];
        for i=1:length(restid)
            id=restid(i);
            %find closest nan to define largest acceptable size of symmetric window
            maxdist = min(abs(nanids-id))-1;
            if maxdist < 2
                tmp(id,:) = inputcol(id,:);
            else
                window = id-maxdist:id+maxdist;
                tmp(id,:) = mean(inputcol(window,:));
            end
        end
        
        outputcol = tmp;
    else
        outputcol = inputcol;
    end
    
    output(:,col) = outputcol;
end

if transposed
    output = output';
end

end