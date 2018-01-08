 function ind=mcnumbers(name,labels)

 %{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be

Extract corresponding labels indices
 
 Use: ind=mcnumbers(d.markerNamee,{'NameOfMarker1','NameOfMarker2','LeftFoot',...});
 
%}
 
ind=zeros(length(labels),1);
for i=1:length(labels)
    ind(i)=find(ismember(name,labels{i}));
end

 end