function PopupProgress(name, varargin)
if nargin == 1
  incrementPopupTicks(name);
elseif (nargin == 2 || nargin == 3) && isfloat(varargin{1})
  makePopupTicks(name, varargin{:})
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makePopupTicks(name, numTotalTicks, updatePeriod)
if nargin < 3
  updatePeriod = 10.0;
end
f = NamedFigure(name);
clf
set(f, 'MenuBar', 'none')

finalStr = sprintf('%lu / %lu', numTotalTicks, numTotalTicks);
width = 1.5 * max(length(finalStr), length(name));

h1 = uicontrol('Style', 'text', 'String', name);
set(h1, 'Units', 'characters')
h1Position = get(h1, 'Position');
h1Position([1 2]) = [0 2];
h1Position(3) = width;
h1Position(4) = 1;
set(h1, 'Position', h1Position)

textStr = sprintf('%lu / %lu', 0, numTotalTicks);
h2 = uicontrol('Style', 'text', 'String', textStr);
set(h2, 'Units', 'characters')
h2Position = get(h2, 'Position');
h2Position([1 2]) = [0 1];
h2Position(3) = width;
h2Position(4) = 1;
set(h2, 'Position', h2Position)
set(h2, 'UserData', 0)  %store numTicks in h2.UserData

h3 = uicontrol('Style', 'text', 'String', '');
set(h3, 'Units', 'characters')
h3Position = get(h3, 'Position');
h3Position([1 2]) = [0 0];
h3Position(3) = width;
h3Position(4) = 1;
set(h3, 'Position', h3Position)
updateStruct.t0 = tic;
updateStruct.updatePeriod = updatePeriod;
updateStruct.nextUpdate = toc(updateStruct.t0) + updatePeriod;
set(h3, 'UserData', updateStruct) %store update schedule in h3.UserData

%store numTotalTicks in f.UserData:
set(f, 'UserData', numTotalTicks)
set(f, 'Units', 'characters')
fPosition = get(f, 'Position');
fPosition(3) = width;
fPosition(4) = 3;
set(f, 'Position', fPosition);
set(f, 'Resize', 'off')
drawnow
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function incrementPopupTicks(name)
f = findobj('Name', name);
children = get(f, 'Children');
h2 = children(2);
numTotalTicks = get(f, 'UserData');
numTicks = get(h2, 'UserData') + 1;
if numTicks == numTotalTicks
  close(f)
  return
end

set(h2, 'UserData', numTicks)
textStr = sprintf('%lu / %lu', numTicks, numTotalTicks);
set(h2, 'String', textStr)

h3 = children(1);
updateStruct = get(h3, 'UserData');
tNow = toc(updateStruct.t0);
if tNow > updateStruct.nextUpdate
  tFinish = (numTotalTicks - numTicks ) * tNow / numTicks;
  timeStr = [getTimeString(tFinish), ' remaining'];
  set(h3, 'String', timeStr)

  updateStruct.nextUpdate = tNow + updateStruct.updatePeriod;
  set(h3, 'UserData', updateStruct)
  
  fPosition = get(f, 'Position');
  minLen = 1.5 * length(timeStr);
  if fPosition(3) < minLen
    fPosition(3) = minLen;
    set(f, 'Position', fPosition)
    for n = 1:length(children)
      pos = get(children(n), 'Position');
      pos(3) = minLen;
      set(children(n), 'Position', pos)
    end
  end
end
drawnow
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timeString = getTimeString(tFinish)
s = round(tFinish);
if s < 60
  timeString = sprintf('%ds', s);
  return
end
m = floor(s/60);
s = s - m * 60;
if m < 60
  timeString = sprintf('%dm%ds', m, s);
  return
end
h = floor(m/60);
m = m - h * 60;
if h < 24
  timeString = sprintf('%dh%dm%ds', h, m, s);
  return
end

d = floor(h/24);
h = h - d * 24;
timeString = sprintf('%gd%dh%dm%ds', d, h, m, s);
return