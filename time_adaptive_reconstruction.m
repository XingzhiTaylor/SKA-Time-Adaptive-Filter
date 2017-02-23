%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This MATLAB function implements the generalized time-adaptive       %
%  reconstruction kernel derived by Yufang Hao and Achim Kempf.            %                                                                                                               
%      The kernel takes in some samples recorded on an unequidistant grid, %
%  and reconstruct the signal on a dense grid. In other words, the function%  
%  interpolates amplitude values that are in between the given samples.    %
%                                                                          %
%      The function takes four input arguments. Here is the introduction   %                                                                                                                                     
%      1.sample: These are the sample amplitudes recorded on the           %    
%  unequidistant grid. The function will interpolate other amplitudes      %
%  between them.                                                           %
%      2.sample_grid: The unequidiatant grid corresponding to the sample   %
%  amplitudes.                                                             %                                      
%      3.longer_sample_grid: This is actually the same grid as             %
%  the sample_grid, but with 500 extra points on both ends.                %
%  When we are evaluate the value of the sine series, we need to add up    %
%  1000 terms centered at a particular t value. However, we cannot find    %
%  all the 1000 terms that center at the beginning and the end of the      %
%  vector. Therefore, we needs these extra points when the center is at    %
%  the beginning and the end of the vector.                                %
%      4.grid_points: This is the new sampling grid that the user wants    %
%  to reconstruct the signal up to.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recons = time_adaptive_reconstruction(sample,sample_grid,longer_sample_grid,grid_points)
    grid_points = round(grid_points,5);                                    
    %Get rid of errors due to the double-precision
    sample_grid = round(sample_grid,5);
    grid_len = length(grid_points);
    recons = zeros(1, length(grid_points));                                
    %Preallocation           
    %% Add more tms on both ends
    longer_tp = t_prime(longer_sample_grid);                               
    %Find the t' values of the sample grid
    %% Reconstruction Part
    for a = 1:grid_len;                                                   
        %For every t in the new grid, find out its amplitude
        t = grid_points(a);                                                
        %Pick a t from the new grid
        k = G(t,sample_grid,longer_sample_grid,longer_tp);                 
        %This is the time-adaptive reconstruction kernel
        recons(a) = dot(sample, k);                                        
        %Multiply each sample value and the reconstruction kernel. Then add them together.
    end
end
%% 
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
