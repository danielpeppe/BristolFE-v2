clear
close all

figure;
xlabel('Time (s)')
ylabel('Magnitude (-)')

hold on
load('g4-s8_x-8000_z0_025.mat');
%Define aperture
aperture_n_els = 8; %number of elements
aperture_vec = 1:aperture_n_els:128;
%Loop over apertures
for i = 1:(numel(aperture_vec) - 1)
    %Get time data for aperture
    aperture = aperture_vec(i):aperture_vec(i + 1);
    aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx,aperture);
    aperture_time_data = sum(exp_data.time_data(:,aperture_data),2);
    %Translate each aperture_time_data to line them up
    large_dsp_i_vec = find(abs(aperture_time_data) > 1);
    translate_time = exp_data.time(large_dsp_i_vec(1));
    %Plot
    plot(exp_data.time - translate_time, aperture_time_data);
end
hold off
