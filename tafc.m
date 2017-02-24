%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   "tafc" stand for time-adaptive filtering & reconstruction. This       % 
% function is the combination of the time-adaptive filter, time-adaptive  %
% reconstruction, and the grid calculation. The tafc function can direct- %
% ly output the signal after filtering that has the same length as the    %
% raw signal. The function needs three vector inputs: raw signal, the     %
% dense sampling grid on which the raw signal is sampled, and the bw vec- %
% tor that holds the bandwidth information of the signal.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtered_signal = tafc(raw_signal,raw_grid,bw)
%% Calculate the sampling grids
    filtered_grid = zeros(1,length(raw_grid));
    %Preallocation. The length of the sparse grid is originally set to be equal
    %to the signal since the real length of the sparse sampling grid is unknown. 
    %Redundant 0 will be discarded at last.
    spacing = raw_grid(2) - raw_grid(1);
    ii = raw_grid(1);
    jj = 1;
    while ii <= raw_grid(end)
        filtered_grid(jj) = ii;
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
        jj = jj + 1;
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
    %If the bandwidth is too low, there may be only one raw_signal point in the
    %filtered_grid. Then the grid should be extended in another way. 
    longer_filtered_grid = horzcat(filtered_grid(1)-500*l:l:filtered_grid(1)-l,filtered_grid,filtered_grid(end)+r:r:filtered_grid(end)+500*r);
    %These are the extra 500 points for the series in the kernel. The grid is
    %extended linearly
    V_t = 2*bw;
%% Filter and down-sample the signal
    raw_grid = round(raw_grid,5);                                                  
    %Get rid of errors due to the double-precision
    filtered_grid = round(filtered_grid,5);
    filtered_signal_short = zeros(1, length(filtered_grid));                             
    %Preallocation
    longer_filtered_tp = t_prime(longer_filtered_grid);
    %Find the t' values of the new grid
    sq = zeros(1,length(raw_grid));                                               
    %Preallocation of sq vector. sq vector stores the value of the sine series. 
    %f = longer_filtered_grid(501:end-500);                                        
    %f is just the new grid.
    %f = round(f,5);                                                 
    %Get rid of errors due to the double-precision
    for ii = 1:length(raw_grid)                                                    
       %For each t_r in the dense sampling grid, evaluate the sine series with respect to t_r
       t_r = raw_grid(ii);                                                                
       sq(ii) = sum_square(t_r,longer_filtered_grid,longer_filtered_tp);           
       %The sum_square function evaluates the value of the sine series
    end
    group = zeros(1,length(filtered_grid)+1);
    %group(1) = length(raw_grid(raw_grid < f(1)));group(ii-1)+ & raw_grid > f(ii-1)
    group(end) = length(raw_grid);
    for ii = 1:length(group)-1
        group(ii) = length(raw_grid(raw_grid <= filtered_grid(ii)));
    end
    sign = ones(1,length(raw_grid));
    for ii = 1:2:length(group)-1
        sign(group(ii)+1:group(ii+1)) = sign(group(ii)+1:group(ii+1)) * (-1);
    end    
    for a = 1:length(filtered_grid)                                          
        %For every t in the new grid, find out its amplitude            
        new_t = filtered_grid(a);                                                  
        %Pick a t from the new grid
        new_tp = longer_filtered_tp(a+500);                                        
        %Pick the corresponding t' value of variable new_t 
        new_t = round(new_t,5);
        b = sqrt(new_tp).*sq./abs(new_t-raw_grid).*V_t*spacing;
        sign = (-1)*(sign);
        if a == 1
            sign(1:group(a)) = sign(1:group(a))*(-1);
        else
            sign(group(a-1)+1:group(a)) = sign(group(a-1)+1:group(a))*(-1);
        end
        %-1 to the power of z gives the sign
        %sign(1) = (-1)*sign(1);
        b = b.*sign;   
        %Assign the sign
        b(isnan(b)|isinf(b)) = 1;
        %Multiply each raw_signal value and the filtering kernel. Then add them together.*0.1  
        filtered_signal_short(a) = dot(raw_signal, b);
        %fprintf('%d\n',a);
        %break;
        %Print out the value of a, indicating which point we are working on. This statement is completely optional.
    end
%% Up-sample the filtered signal
    longer_tp = t_prime(longer_filtered_grid); 
    filtered_signal = zeros(1,length(raw_grid));
    for a = 1:length(raw_grid);                                                   
        %For every t in the new grid, find out its amplitude
        t = raw_grid(a);                                                
        %Pick a t from the new grid
        k = G(t,filtered_grid,longer_filtered_grid,longer_tp);                 
        %This is the time-adaptive filtered_signaltruction kernel
        filtered_signal(a) = dot(filtered_signal_short, k);                                        
        %Multiply each sample value and the filtered_signaltruction kernel. Then add them together.
    end
end
function a = G(t,tn,longer_tm,longer_tp)                                   
    %This is the reconstruction kernel
    sq = sum_square(t,longer_tm,longer_tp);                                                    
    a = (sqrt(longer_tp(501:end-500))./abs(t-tn))*sq;   
    a(isnan(a)|isinf(a)) = 1;
    sign = zeros(1,length(tn));
    len = length(tn(tn<t))-1;
    sign(1:len+1) = len:-1:0;
    sign(len+2:end) = 0:length(tn)-2-len;
    sign_1 = (-1).^(sign);
    a = a.*sign_1;
end
%% 
function sq = sum_square(t,longer_tm,longer_tp)                            
%Calculate the sine series
   center = length(longer_tm(longer_tm<t));                                
   %For the series, we need to find 1000 tms that center at t_hat. We first find the index of the center
   tm = longer_tm(center-499:center+500);                                  
   %Pick those tms according to the center
   tm_p = longer_tp(center-499:center+500);                                
   %Pick the corresponding tm'
   s = 1./(t-tm); s_1 = s.*tm_p;                                           
   %Store terms in the series in two row vectors. They are 1/(t_hat-tm) and tm'/(t_hat-tm)
   sq = (dot(s,s_1)).^(-1/2);                                              
   %Sum them up using dot product
end
%% 
function t_p = t_prime(t)
    t_vector = zeros(1,length(t)+2);                                       
    %Calculate the t' values
    t_vector(1) = t(1); t_vector(end) = t(end);
    t_vector(2:end-1) = t;
    difference = t_vector(3:end)-t_vector(1:end-2);                        
    %Take the spacing between sample points as t'
    t_p = difference/2;
end
