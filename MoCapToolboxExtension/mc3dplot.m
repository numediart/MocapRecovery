function mc3dplot(d,p,myfighandle,speed)

global isplaying;
global gframeid;
global gspeed;
gspeed = 10;
gframeid = 1;
isplaying = false;
%{
Created by Mickaël Tits, numediart Institute, University of Mons, Belgium
21/12/2017
Contact: mickaeltits@gmail.com or mickael.tits@umons.ac.be
 
3D display of a MoCap structure from the MoCap Toolbox

 inputs:
 - d: MoCap structure from the MoCap Toolbox
 - p: MoCap animation parameters (see mcinitanimpar)
 - myfighandle: handle to a figure (optional)
%}

if nargin < 2 || isempty(p)
   p = mcinitanimpar(); 
end
if nargin < 3
    myfighandle = figure;
end

if nargin == 4
   gspeed = speed; 
end

try
    fig = myfighandle;
catch
    fig = figure;
end

clf;

allx = d.data(1:3:end,:);
ally = d.data(2:3:end,:);
allz = d.data(3:3:end,:);

p.minx = nanmin(nanmin(allx));
p.miny = nanmin(nanmin(ally));
p.minz = nanmin(nanmin(allz));
p.maxx = nanmax(nanmax(allx));
p.maxy = nanmax(nanmax(ally));
p.maxz = nanmax(nanmax(allz));

set(fig,'Position',[50 50 p.scrsize(1) p.scrsize(2)]);

fig.CloseRequestFcn = @myclose_req;

scatter3(1,1,1);%dummy plot to initialize view orientation
plotframe(d,p,gframeid);

hold on;
uislider = uicontrol('Parent',fig,'Style','slider','Position',[0 20 p.scrsize(1) 20],...
    'value',1, 'min',1, 'max',d.nFrames);
uislider.Callback = @(es,ed) plotframe(d,p,es.Value);

uiplaystop =  uicontrol('Style', 'togglebutton', 'String', 'Play__Pause',...
    'Position', [20 40 100 20]);
uiplaystop.Value = isplaying;
uiplaystop.Callback = @(es,ed) playloop(d,p,es.Value);

speedslider = uicontrol('Parent',fig,'Style','slider','Position',[130 40 200 20],...
    'value',gspeed, 'min',1, 'max',200);
speedslider.Callback = @(es,ed) setspeed(es.Value);


%uislider.HandleVisibility = 'off';
hold off;
playloop(d,p,isplaying);
end

%%

function plotframe(d,p,newframeid)

global isplaying;
global gframeid;
gframeid = newframeid;
if gframeid > d.nFrames
    gframeid = 1;
    newframeid = 1;
end
[az,el] = view();

newframeid = round(newframeid);
if newframeid < 1
    newframeid = 1;
elseif newframeid > d.nFrames
    newframeid = d.nFrames;
end

frame = d.data(newframeid,:);
x = frame(1:3:end);
y = frame(2:3:end);
z = frame(3:3:end);
hold off;
if isempty(p.markercolors)
    scatter3(x,y,z,'b');
else
    scatter3(x,y,z,36,p.markercolors);
end
hold on;
try
    if isempty(p.conncolors)
        for i = 1:length(p.conn)
            plot3(x(p.conn(i,:)),y(p.conn(i,:)),z(p.conn(i,:)),'b');
        end
    else
        for i = 1:length(p.conn)
            plot3(x(p.conn(i,:)),y(p.conn(i,:)),z(p.conn(i,:)),'Color',p.conncolors(i,:));
        end
    end
catch
    warning('Connections data does not correspond to markers data. Check your parameter p.conn, or simply remove it.');
end


axis([p.minx p.maxx p.miny p.maxy p.minz p.maxz]);
view(az,el);
try
    drawnow;
catch
    isplaying = false;
end
end

function playloop(d,p,playing)

global isplaying;
global gframeid;
global gspeed;
isplaying = playing;

while isplaying
    gframeid = gframeid+round(gspeed);
    if gframeid > d.nFrames
        gframeid = 1;
    end
    try
        fig = gcf;
    catch
        isplaying = false;
        return;
        %Figure probably closed without stopping the player
    end
    try
        fig.Children(length(fig.Children)-1).Value = gframeid;
        plotframe(d,p,gframeid)
    catch
        %warning('UI slider could not be updated. Figure probably closed without stopping the player.');
        isplaying = false;
    end
end

end

function setspeed(speed)

global gspeed;
gspeed = speed;
end

function myclose_req(src,callbackdata)

global isplaying;
isplaying = false;
closereq;

end

