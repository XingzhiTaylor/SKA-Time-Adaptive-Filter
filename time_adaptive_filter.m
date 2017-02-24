%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This MATLAB function implements the generalized time-adaptive       %
%  filter derived by Yufang Hao and Achim Kempf.                           %                                                                                                               
%      The filter takes in a raw signal as amplitudes recorded on a dense  % 
%  sampling grid. It also takes the information of the new sampling grid   %    
%  that the user wants to filter the signal down to. The output is the     %
%  filtered signal given by amplitudes on the new sampling grid.           %                                                                                                                           
%                                                                          %
%      The function takes five input arguments. Here is the introduction   %                                                                                                                                     
%      1.sample: This is the raw signal before filtering. The vector       %
%  constains all the amplitudes on the dense grid                          %
%      2.raw_grid: This is the old dense sampling grid called {t_n}        % 
%      3.filtered_grid: This is the new sampling grid that the             %
%  user wants to filter the signal down to. This grid is called {t_r}.     %
%      4.longer_filtered_grid: This is actually the same grid as           %
%  the filtered_grid, but with 500 extra points on both ends.              %
%  When we are evaluate the value of the sine series, we need to add up    %
%  1000 terms centered at a particular t value. However, we cannot find    %
%  all the 1000 terms that center at the beginning and the end of the      %
%  vector. Therefore, we needs these extra points when the center is at    %
%  the beginning and the end of the vector.                                %
%      5.V_t: This vector contains the bandwidth information at each       %
%  point on the dense grid. According to the formula, we have to multiply  %
%  each sample amplitude (those on the dense grid) with the two times the  %
%  bandwidth at that point.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtered_signal = time_adaptive_filter(sample,raw_grid,filtered_grid,longer_filtered_grid,V_t)
    raw_grid = round(raw_grid,5);                                                  
    %Get rid of errors due to the double-precision
    filtered_grid = round(filtered_grid,5);
    spacing = raw_grid(2) - raw_grid(1);
    filtered_signal = zeros(1, length(filtered_grid));                             
    %Preallocation
    longer_filtered_tp = t_prime(longer_filtered_grid);
    %Find the t' values of the new grid
    sq = zeros(1,length(raw_grid));                                               
    %Preallocation of sq vector. sq vector stores the value of the sine series. 
    for ii = 1:length(raw_grid)                                                    
       %For each t_r in the dense sampling grid, evaluate the sine series with respect to t_r
       t_r = raw_grid(ii);                                                                
       sq(ii) = sum_square(t_r,longer_filtered_grid,longer_filtered_tp);           
       %The sum_square function evaluates the value of the sine series
    end
    %The sign of a sample point in {t_r} grid is determined by the number of points in {t_r} grid that lie between the 
    %t_r point being calculated and t_n point being read. The t_n points that lie between the same two t_r points will
    %give the same signs for the t_r point being calculated. Therefore, dividing {t_n} into groups makes the sign calculation
    %faster.
    group = zeros(1,length(filtered_grid)+1);
    group(1) = length(raw_grid(raw_grid < f(1)));
    group(end) = length(raw_grid);
    %The vector "group" holds the information of starting and ending indeices of each {t_n} group
    for ii = 2:length(group)-1
        group(ii) = group(ii-1)+length(raw_grid(raw_grid < f(ii) & raw_grid > f(ii-1)));
    end
    %Then, create a vector that holds the signs for {t_n}. 
    sign = ones(1,length(raw_grid));
    for ii = 1:2:length(group)-1
        sign(group(ii)+1:group(ii+1)) = sign(group(ii)+1:group(ii+1)) * (-1);
    end
    %Initialize the sign vector
    for a = 1:length(filtered_grid);                                               
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
            %Calculate the sign vector in terms of the first t_r
        else
            sign(group(a-1)+1:group(a)) = sign(group(a-1)+1:group(a))*(-1);
            %The signs do not need to be repeatedly calculated. We just need to flip them.
        end
        %-1 to the power of z gives the sign   
        b = b.*sign;   
        %Assign the sign
        b(isnan(b)|isinf(b)) = 1;
        %Multiply each sample value and the filtering kernel. Then add them together.*0.1  
        filtered_signal(a) = dot(sample, b);
        fprintf('%d\n',a);
        %Print out the value of a, indicating which point we are working on. This statement is completely optional.
    end
end
%%
function sq = sum_square(t,longer_tm,longer_tp)                                     
%Calculate the sine series
   center = length(longer_tm(longer_tm<t));                                         
   %For the series, we need to find 1000 tms that center at t_hat. We first find the index of the center
   tm = longer_tm(center-499:center+499);                                           
   %Pick those tms according to the center
   tm_p = longer_tp(center-499:center+499);                                         
   %Pick the corresponding tm'
   s = 1./(t-tm); s_1 = s.*tm_p;                                                   
   %Store terms in the series in two row vectors. They are 1/(t_hat-tm) and tm'/(t_hat-tm)
   sq = (s*s_1').^(-1/2);                                                           
   %Sum them up using dot product
end
%% 
function t_p = t_prime(t)                                                           
    %Calculate the t' values
    t_p = zeros(1,length(t));                 
    t_p(1:end-1) = diff(t);                                                        
    %Take the spacing between sample points as t'
    t_p(end) = t_p(end-1);                               
end 
