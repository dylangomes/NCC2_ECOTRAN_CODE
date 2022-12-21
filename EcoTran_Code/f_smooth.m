function smeanx = f_smooth(x, p)
% preform mean smoothing of data
% by Jim Ruzicka
% x is data array
% smoothing done along rows for each column of data
% p is number of points before and after the target point for smoothing
% revision date: 3-12-2014

% determine size of data array
[numrows, numclms] = size(x);

for loop1 = 1:numrows
    
    windowstart = loop1 - p;
    windowend   = loop1 + p;
    
    if windowstart < 1 
        windowstart = 1; 
    end % truncate window at end points
    
    if windowend > numrows
        windowend = numrows;
    end
    
    window = x(windowstart:windowend, :);
    
    % look to see if window is all NaN
    look = find(~isnan(window));
    
    if isempty(look)
        smoothwindow = NaN;
    else
        smoothwindow = mean(window); % change this to nanmean as soon as stat toolbox available
    end
    
    smeanx(loop1, 1:numclms) = smoothwindow;
    
end % end loop1

% end m-file***************************************************************