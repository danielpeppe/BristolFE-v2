function fn_plot_porosity_correlations(op_save, fieldname, save_figure, dir_name)

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
ylabel('Attenuation [% of front-wall signal]')
xlabel(string(fieldname))

%% Save figure
if save_figure
    print(figure1,'-dpng',['-r','1000'], [dir_name '/' fieldname '_plot.png'])
end


end