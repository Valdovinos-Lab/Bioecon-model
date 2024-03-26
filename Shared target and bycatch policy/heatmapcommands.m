%%Commands to generate heatmaps for Figures 3 and 4
% 1) Median fraction of Target surviving (yr 5 of fishing/penultimate fishing free year) 
% 2) Median fraction of Bycatch surviving (yr 5 of fishing/penultimate fishing free year)
% 3) Mean target extinctions (since fishing started)
% 4) Mean bycatch extinctions (since fishing started)
% 5) Mean secondary extinctions (since fishing started)
% 6) Mean target yield (year 5 of fishing)
% 7) Mean bycatch yield (year 5 of fishing)
% 8) Mean profits (year 5 of fishing)
% 9) average fishing year (year 5 of fishing)

xvalues =  [0.15 0.3 0.45 0.6];
yvalues = [0.6 0.45 0.3 0.15 0];
data = zeros(5,4);

%aj = policydata{7};  %random bycatch
%aj2 = policydata{6};

aj = policydata_min{7};  %vulnerable bycatch
aj2 = policydata_min{6};
for i=1:4
    for j=1:5
        data(j,i) = aj(i,6-j)/aj2(i,6-j); %calculating bycatch per unit target Fig 4F
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
title('Average bycatch/target (year 5)');
xlabel('esc_{target}');
ylabel('esc_{bycatch}');
