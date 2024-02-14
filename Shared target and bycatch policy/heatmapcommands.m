xvalues =  [0.15 0.3 0.45 0.6];
yvalues = [0.6 0.45 0.3 0.15 0];
data = zeros(5,4);
%aj = policydata_vulnerablebycatch{5}; %3 - target, 4-bycatch, 5-sec
aj = policydata_min{4}; %1 - median target, 2 - median bycatch, 8 - profits, 7/6 is bycatch/target yield
for i=1:4
    for j=1:5
        data(j,i) = aj(i,6-j);
    end
end
zvalues = data;
x_finer = linspace(0.15, 0.6, 50);
y_finer = linspace(0.6,0,50);
[x_finer_grid, y_finer_grid] = meshgrid(x_finer, y_finer);
z_finer_matrix = interp2(xvalues, yvalues, zvalues, x_finer_grid, y_finer_grid, 'spline');
%figure
subplot(2,3,5)

h = heatmap(x_finer, y_finer, z_finer_matrix, 'GridVisible', 'off');
newXLabels = cell(1,50);
newYLabels = cell(1,50);
for i=1:50
    newXLabels{i} = '';
    newYLabels{i} = '';
end
newXLabels{1} = 0.15;
newXLabels{17} = 0.3;
newXLabels{35} = 0.45;
newXLabels{50} = 0.6;
newYLabels{1} = 0.6;
newYLabels{13} = 0.45;
newYLabels{25} = 0.3;
newYLabels{38} = 0.15;
newYLabels{50} = 0;
% Update the X-axis and Y-axis labels
h.XDisplayLabels = newXLabels;
h.YDisplayLabels = newYLabels;
title('');
xlabel('esc_{target}');
ylabel('Average profits (year 5)');
