function missing_markers=mcviewmissing(d,display)
 
 %{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Extract and visualize missing data in a MoCap structure
 
%}

[mf, mm, mgrid] = mcmissing(d);
missing_markers=d.markerName(find(mf));

if nargin < 2
   display = 1; 
end

if display
    disp('missing markers :');
    missing_markers
    figure;
    subplot(3,1,1), bar(mf), xlabel('Marker'), ylabel('Num. of Missing frames')
    subplot(3,1,2), bar(mm), xlabel('Frame'), ylabel('Num. of Missing markers')
    subplot(3,1,3), imagesc(-mgrid'), colormap gray, xlabel('Frame'),...
    ylabel('Marker')
end

 end