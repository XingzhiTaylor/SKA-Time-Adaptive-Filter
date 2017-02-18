function signal = time_varying_signal_generator(grid,longer_grid)
    %The sample signal generator takes a sampling grid as a input
    %since the sampling grid contains all the information of its 
    %corresponding time-varying bandwidth
    longer_tp = t_prime(longer_grid);
    %% Construct the signal
    signal = zeros(1,length(grid));
    factor = randn(1,length(grid));                                              
    %The amplitude of signal has a Gaussian distribution
    for a = 1:length(grid)
       t = grid(a);
       k = G(t,grid,longer_grid,longer_tp);                                          
       %The signal consists of many generalized sinc functions that center at each tn
       signal(a) = dot(factor,k);
    end
end
%%
function a = G(t,tn,longer_grid,longer_tp)
    sq = sum_square(t,longer_grid,longer_tp);                        
    a = (sqrt(longer_tp(501:end-500))./abs(t-tn))*sq;   
    a(isnan(a)|isinf(a)) = 1;
    sign = zeros(1,length(tn));
    len = length(tn(tn<t))-1;
    sign(1:len+1) = len:-1:0;
    sign(len+2:end) = 0:length(tn)-2-len;
    sign = (-1).^(sign);
    a = a.*sign;
end
%% 
function sq = sum_square(t,longer_grid,longer_tp)       
   %Calculate the series
   center = length(longer_grid(longer_grid<t));   
   tm = longer_grid(center-499:center+500);
   tm_p = longer_tp(center-499:center+500);
   s = 1./(t-tm); s_1 = s.*tm_p;
   sq = sqrt(1./(s*s_1'));                            
   %Sum them up
end
%% 
function t_p = t_prime(t)
    t_vector = zeros(1,length(t)+2);
    t_vector(1) = t(1); t_vector(end) = t(end);
    t_vector(2:end-1) = t;
    difference = t_vector(3:end)-t_vector(1:end-2);
    t_p = difference/2;
end