function fn_plot_results(op, op_output, res, steps, exp_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Plot Final Results
fig = figure;
FntN = 'Cambria';
FntS = 16;

%Get scale factor to scale exp data (by first responses)
if op_output.plot_exp_data
    res_sum_dsps = sum(res{1,1}{1}.dsps);
    %exp aperture
    aperture = 1:op.aperture_n_els;
    aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
    aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
    %Scale exp response to match sim response
    scale_exp = max(abs(aperture_dsp_data))/max(abs(res_sum_dsps));
end

%Plot simulated data
if op_output.plot_sim_data
    %Iterate over sim results
    for i = 1:length(res)
        %sum over all the nodes where displacments were recorded
        res_sum_dsps = sum(res{1,i}{1}.dsps);
    
        %Get params for legend when running multiple sims
        if ~isempty(op.params)
            if iscell(op.params)
                params_i = op.params{i};
                n = length(params_i);
                if n == 2
                    for ii = 1:(n - 1)
                        params_i = strcat(string(params_i(1)),string(params_i(2)));
                    end
                end
                str = params_i;
            else
                str = string(op.params(i));
            end
        else
            str = 'Simulated response';
        end
    
        %Plot final results
        plot(steps{1,i}{1}.load.time, res_sum_dsps,'DisplayName',str, 'Color', 'k', 'LineStyle','-', 'LineWidth', 1);
        hold on
    end
end

% Plot experimental data on top
if op_output.plot_exp_data
    hold on
    %Plot
    translate_time = 1.1955e-05; %Start of exp response
    plot(exp_data.time - translate_time, aperture_dsp_data/scale_exp, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','Experimental response','LineWidth', 1);
    hold off
end

xlabel('Time [s]')
ylabel('Displacement [-]')

xcorr = 2; ycorr = 2;
w = 30; h = 15;
osx = 2; osy = 2;
% osx = 2; osy = 3.5;
set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS)
legend('Location','north'); legend boxoff

end

