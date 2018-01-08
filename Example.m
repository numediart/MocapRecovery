
%This example requires the MoCap Toolbox available at :
% https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox
% or the direct link to avoid giving your e-mail :
% https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox/@@resolveuid/195072f40a274f57a893387de8ecd695

% Download the toolbox and add the folder 'mocaptoolbox' to the same folder as this script
% Add the data you want to test to the same folder as this script

%% Initialize and load data
clear all;
close all;
clc;

%Do not forget to download and add the mocap toolbox folder in this repository
addpath('mocaptoolbox');
addpath('mocaptoolbox/private');
addpath('MoCapToolboxExtension'); %some extensions to the MoCap Toolbox

connectionFile = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can change these parameters

%input files (marker data and connections)
%dataFile = '85_02.c3d';
connectionFile = [];

 dataFile = '01_karate.c3d';
 connectionFile = '01_karate.txt';

%dataFile = '03_rolling.c3d';
%connectionFile = '03_rolling.txt';

%dataFile = '02_dancingandfalling.c3d';
%connectionFile = '02_dancingandfalling.txt';



use_manual_connection = true;

if isempty(connectionFile)
    use_manual_connection = false;
end

% if you don't have any connection file, put use_manual_connection to
% false. Automatic bones (connections) will be designed instead

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data (requires the MoCap toolbox)
original = mcread(['data/' dataFile]);

%Remove invalid markers in hdm files
markernames = original.markerName;
ids = find(~ismember(markernames,{'*0','*1','*2'}));
original = mcgetmarker(original,ids);

%interpolate some unexpected really small gaps in some hdm files
if strcmp(lower(original.filename(1:3)),'hdm')
    numberofnans = sum(sum(isnan(original.data)));
    fprintf('number of nans in original file: %d\n',numberofnans);
    if numberofnans > 0
        figure(1);imagesc(isnan(original.data));
    end
    original.data = naninterp(original.data,'pchip');
    disp('salut');
end

%% Display original data

%MoCap display parameters
p = mcinitanimpar;
%manual connections (from file)
if use_manual_connection
    p=mccreateconnmatrix(['data/' connectionFile],p);
    fprintf('Using manual connection matrix...\n');
else
    %or automatic connections
    %p.conn = mcautomaticbones(original,30,400,20);
    p.conn = mcautomaticbones2(original);    
    fprintf('Using automatic connection matrix...\n');
end
%figure window size
p.scrsize = [1900/2 800];

%Display original data
myfighandle = figure(1);
mc3dplot(original,p,myfighandle);
title('Original data');

%%% If you want to recover from an incomplete original sequence:

incomplete = original;

%% (Optional) you can use this cell to simulate gaps on your data:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can change these parameters
nmissing = 5; %number of missing markers
gapsec = 2; %duration of gaps
multiple_gaps = true; %if true, multiple gaps in a marker trajectory (random number)
max_number_of_multiple_gaps = 10;% maximum possible number of multiple gaps in a marker trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if any(any(isnan(original.data)))
    
    %Original sequence is already incomplete
    warning(['Note that there are already missing data in your file.' ...
        ' Dont use this cell if you want to recover missing data from the original file.']);
end

%Simulate gaps

incomplete = original;
indices = zeros(size(incomplete.data));
nframes=incomplete.nFrames;
nmarkers = incomplete.nMarkers;
fps = incomplete.freq;
j=randperm(nmarkers);
j=j(1:nmissing);
gapframe=round(gapsec*fps);
%%%%one gap per marker
gaps=ones(nmarkers,1);

if multiple_gaps
    %%%%several gaps per marker
    gaps=round(rand(nmarkers,1)*max_number_of_multiple_gaps);
    gaps(gaps<1)=1;
end

indicesperm=[];
%gaps at the same time
id1 = 1+floor(rand(1)*(nframes - gapsec*fps));
id2 = min(nframes,id1+gapframe);
for m=j
    indicesperm{m} = zeros(incomplete.nFrames,1);
    for g=1:gaps(m)
        if multiple_gaps
            %gaps at different times
            id1 = 1+round(rand(1)*nframes*4/5.1);
            id2 = min(nframes,id1+gapframe);
        end
        indicesperm{m}(id1:id2)=true;
        indices(id1:id2,(3*m-2):(3*m))=true;%indexes of fake missing frames
    end
    indicesperm{m}=logical(indicesperm{m});
end
indices = logical(indices);
incomplete.data(indices)=nan;

figure(2);imagesc(indices);
title('Simulated gaps');
xlabel('Markers (interleaved x,y,z)');
ylabel('Frames');




%Display original data
myfighandle = figure(1);
mc3dplot(incomplete,p,myfighandle);
title('Fake data');

%% Recovery examples
%% Gloersen et al. 2016 only

%Check all possible options in the mcrecovery function
%options.saveastsv = 1;%save recovered sequence as tsv file
%options.recursivefilling = 0;%use or not recusrivefilling
%choose which method to use
options.method1 = 0;%Local interpolation
options.method2 = 0;% Local polynomial regression
options.method3 = 0;%Local GRNN
options.method4 = 0;%Global weighted linear regression
options.method5 = 1;%Gloersen et al. 2016
options.advancedordering = 1;
options.spaceconstraint = 0;%use or not spaceconstraint
options.timeconstraint = 0;%use or not timeconstraint
options.filtering = 0;%use or not timeconstraint
options.quiet = 0;%avoid console output
options.presenceMin = 30;%threshold (in % of available frames) under which discard some markers
% Gloersen would fail without discarding too bad markers

tic;
recovered = mcrecovery(incomplete,options);
toc;

myfighandle = figure(3);
mc3dplot(recovered,p,myfighandle);
title('PCA (Gloersen et al.)');

%% Gloersen et al. 2016 with input filtering and soft constraints

%clc;

%Check all possible options in the mcrecovery function
%options.saveastsv = 1;%save recovered sequence as tsv file
%options.recursivefilling = 0;%use or not recusrivefilling
%choose which method to use
options.method1 = 0;%Local interpolation
options.method2 = 0;% Local polynomial regression
options.method3 = 0;%Local GRNN
options.method4 = 0;%Global weighted linear regression
options.method5 = 1;%Gloersen et al. 2016
options.advancedordering = 1;
options.spaceconstraint = 1;%use or not spaceconstraint
options.timeconstraint = 1;%use or not timeconstraint
options.filtering = 1;%use or not timeconstraint
options.quiet = 1;%avoid console output
options.presenceMin = 30;%threshold (in % of available frames) under which discard some markers

tic;
recovered = mcrecovery(incomplete,options);
toc;

myfighandle = figure(4);
mc3dplot(recovered,p,myfighandle);
title('PCA (Gloersen et al.) with filtering and constraints');

%% Local GRNN

%Different combinations may be tried to extract the best results

%options.saveastsv = 0;%save recovered sequence as tsv file
%options.recursivefilling = 0;%use or not recusrivefilling
%choose which method to use
options.method1 = 0;%Local interpolation
options.method2 = 0;% Local polynomial regression
options.method3 = 1;%Local GRNN
options.method4 = 0;%Global weighted linear regression
options.method5 = 0;%Gloersen et al. 2016
options.advancedordering = 1;
options.spaceconstraint = 1;%use or not spaceconstraint
options.timeconstraint = 1;%use or not timeconstraint
options.filtering = 1;%use or not timeconstraint
options.quiet = 1;%avoid console output
options.presenceMin = 30;%threshold (in % of available frames) under which discard some markers

tic;
recovered = mcrecovery(incomplete,options);
toc;

myfighandle = figure(5);
mc3dplot(recovered,p,myfighandle);
title('Local GRNN');

%% PMA of different individual models

%Different combinations may be tried to extract the best results

%options.saveastsv = 0;%save recovered sequence as tsv file
%options.recursivefilling = 0;%use or not recusrivefilling
%choose which method to use
options.method1 = 1;%Local interpolation
options.method2 = 1;% Local polynomial regression
options.method3 = 1;%Local GRNN
options.method4 = 1;%Global weighted linear regression
options.method5 = 1;%Gloersen et al. 2016
options.advancedordering = 1;
options.spaceconstraint = 1;%use or not spaceconstraint
options.timeconstraint = 1;%use or not timeconstraint
options.filtering = 1;%use or not timeconstraint
options.quiet = 1;%avoid console output
options.presenceMin = 30;%threshold (in % of available frames) under which discard some markers

tic;
recovered = mcrecovery(incomplete,options);
toc;

myfighandle = figure(6);
mc3dplot(recovered,p,myfighandle);
title('PMA of different individual models');

%%
close all
