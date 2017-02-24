%There should be two row vectors prepared before running this script: 
%     1.The dense sampling grid "raw_grid"
%     2.The bandwidth vector "bw" that stores the bandwidth value for all the
%     points in raw_grid. The bandwidth vector should have the same length
%     as raw_grid and the default unit for bandwidth is Hz.
calculate_grid
%This is the script for grid calculation that generates filtered_grid,
%longer_filtered_grid, and V_t
white_noise = zeros(1,length(raw_grid));
%Pre-allocation of the white_noise vector
for ii = raw_grid
    components = randn*sinc((raw_grid)-ii);
    %The sinc components center at different positions in the raw_grid with
    %Gaussian distributed amplitudes.
    white_noise = white_noise + components;
end
filtered_white_noise_short = time_adaptive_filter(white_noise,raw_grid,filtered_grid,longer_filtered_grid,V_t);
%This is the down-sampled signal. This signal has a small data volumn 
filtered_white_noise_long = time_adaptive_filter(filtered_white_noise_short,raw_grid,filtered_grid,longer_filtered_grid);
%This is the time-adaptively reconstructed signal. This signal has the same length as the white noise and is the final output. 
