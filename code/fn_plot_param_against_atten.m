function fn_plot_param_against_atten(op_save, fieldname, batch_path,  save_figure)

%Plot Final Results
if save_figure
    figure1 = figure('Visible', 'off');
else
    figure1 = figure;
end

x = zeros(size(op_save));
y = zeros(size(op_save));
for i = 1:numel(op_save)
    x(i) = op_save{i}.(fieldname);
    y(i) = op_save{i}.attenuation;
end
[x_sorted, x_i] = sort(x);
y_sorted = y(x_i);

plot(x_sorted, y_sorted)
ylabel('Attenuation coefficient, $\alpha$ [dB/mm]', 'Interpreter', 'latex')
xlabel(string(fieldname))

%% Save figure
if save_figure
    porosity_plot_path = batch_path + "/" + fieldname + "_plot.png";
    print(figure1,'-dpng',['-r','1000'], porosity_plot_path)
end

end