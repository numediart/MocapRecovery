
 function mcwritetsv2(d, path)
 
%{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Modified (bug correction) from the MoCap Toolbox original function mcwritetsv:
Added tests on the existence on optional fields of the input mocap structure d
 
%}



% Saves mocap structure as tsv file.
%
% syntax
% mcwritetsv(d, path)
%
% input parameters
% d: MoCap data structure
% path: path to save the tsv file (optional). If no path is given, file is saved to current directory
%
% output
% tsv file, saved in the current or in the specified directory
%
% examples
% mcwritetsv(d)
% mcwritetsv(d, 'folder')
% mcwritetsv(d, '/path/folder') %(Mac)
%
% comments
% Thanks to Roberto Rovegno, University of Genoa, Italy for contributing to this code
%
% see also
% mcread
% 
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland
%

if ~strcmp(d.type, 'MoCap data')
    disp([10, 'Input is not a MoCap data structure. No file written.' 10])
    return
end

%create file name
[~, name, ~]=fileparts(d.filename);
% name=strtok(name,'.')
name=strcat(name,'.tsv');

if nargin==1
    path=cd;
end

if ~strcmp(path(end),'/')
    path=strcat(path,'/');
end

name = strcat(path,name);

fid = fopen(name, 'w');
fprintf(fid, '%s\t', 'NO_OF_FRAMES');
fprintf(fid, '%u\n', d.nFrames);
fprintf(fid, '%s\t', 'NO_OF_CAMERAS');

test = exist('d.nCameras','var');
if test
    fprintf(fid, '%u\n', d.nCameras);
else
    fprintf(fid, ' --\n'); 
end

fprintf(fid, '%s\t', 'NO_OF_MARKERS');
fprintf(fid, '%u\n', d.nMarkers);
fprintf(fid, '%s\t', 'FREQUENCY');
fprintf(fid, '%u\n', d.freq);
fprintf(fid, '%s\t', 'NO_OF_ANALOG');
fprintf(fid, '%u\n', 0);%d.nAnalog
fprintf(fid, '%s\t', 'ANALOG_FREQUENCY');
fprintf(fid, '%u\n', 0);%d.anaFreq

fprintf(fid, '%s\t', 'DESCRIPTION');
test = exist('d.other.descr','var');
if test
    if strcmp('DESCRIPTION', d.other.descr(1:11))
        fprintf(fid, '%s\n', d.other.descr(12:end));
    else
        fprintf(fid, '%s\n', d.other.descr);
    end
else
    fprintf(fid, ' --\n'); 
end

fprintf(fid, '%s\t', 'TIME_STAMP');
test = exist('d.other.timeStamp','var');
if test    
    if strcmp('TIME_STAMP', d.other.timeStamp(1:10))
        fprintf(fid, '%s\n', d.other.timeStamp(12:end));
    else
        fprintf(fid, '%s\n', d.other.timeStamp);
    end
else
    fprintf(fid, ' --\n'); 
end

fprintf(fid, '%s\t', 'DATA_INCLUDED');
test = exist('d.other.dataIncluded','var');
if test 
    fprintf(fid, '%s\n', d.other.dataIncluded); 
else
    fprintf(fid, ' 3D\n'); %by default, we write 3D data
end

fprintf(fid, '%s\t', 'MARKER_NAMES');
fprintf(fid, '%s\t', d.markerName{:});
fprintf(fid, '%f\n', []);

for k=1:size(d.data,1);
    fprintf(fid, '%f\t', d.data(k,:));
    fprintf(fid, '%f\n', []);
end

fclose(fid);

if exist(name,'file') == 2
    disp([10, name, ' written.', 10])
else
    disp([10, name, ' not written!', 10])
end
   
 end