filtered_grid = zeros(1,length(raw_grid));
%Preallocation. The length of the sparse grid is originally set to be equal
%to the signal since the real length of the sparse sampling grid is unknown. 
%Redundant 0 will be discarded at last.
ii = raw_grid(1);
while ii <= raw_grid(end)
    filtered_grid(fix(ii)) = ii;
    %Store the value of the sampling point
    interp_bw = interp1(raw_grid,bw,ii);
    %Using the built-in linear interpolation function to calculate
    %bandwidth
    space = 1/(2*interp_bw);
    %The spacing should be 1/2*Omega. 
    temp = ii + space;
    if temp == fix(temp)
        ii = temp + 0.01;
    else
        ii = temp;
    end
    %In MATLAB, some number such as 9.999 and 10.001 will be automatically
    %rounded to an integer. The automatic rounding will cause large error.
    %So if an automatic rounding is detected, add 0.01 to change the
    %variable type to a floating point number.
end
filtered_grid = filtered_grid(filtered_grid~=0);
%Discard all 0 values
if length(filtered_grid) == 1
    l = space;
    r = space;
else
    l = abs(filtered_grid(2)-filtered_grid(1));
    r = abs(filtered_grid(end)-filtered_grid(end-1));
end
%If the bandwidth is too low, there may be only one sample point in the
%filtered_grid. Then the grid should be extended in another way. 
longer_filtered_grid = horzcat(filtered_grid(1)-500*l:l:filtered_grid(1)-l,filtered_grid,filtered_grid(end)+r:r:filtered_grid(end)+500*r);
%These are the extra 500 points for the series in the kernel. The grid is
%extended linearly
V_t = 2*bw;