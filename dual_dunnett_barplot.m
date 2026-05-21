function dual_dunnett_barplot(data,genotypes,concentrations,cmap,upperYLabel,defaultsettings)
% Setup Vars and Get Means of Input Data for Bars
num_gt = length(genotypes);
num_conc = length(concentrations);
if nargin<6
    defaultsettings = [
    {[3,1]}, % tiledlayout size 
    {[1,1]}, % upper bar plot tile size
    {[2,1]}, % lower annotations tile size
    {abyss(num_conc)}, % annotations colormap
    ];
end
tl = defaultsettings{1}; % tiledlayout size
uptl = defaultsettings{2}; % upper tiledlayout size
lwrtl = defaultsettings{3}; % lower tiledlayout size
cmap_annotations = defaultsettings{4}; % annotations colormap

numcomparisons = size(data,2)*size(data,3)-size(data,3) + (size(data,3)-1)*size(data,2);


% Input Validation
if (length(genotypes)~=size(data,3)) error("Mismatch between size of genotypes array and dim 3 of input data..."); end
if (length(concentrations)~=size(data,2)) error("Mismatch between size of concentrations array and dim 2 of input data..."); end
if (size(cmap,1)~=size(data,2)) error("Mismatch between colormap and num of concentrations in input data..."); end


plot_means = squeeze(mean(data, 1, "omitmissing"));
plot_sems = squeeze(std(data, 0, 1,"omitmissing")) ./ squeeze(sqrt(sum(~isnan(data),1)));

% Position Axes
tiledlayout(tl(1),tl(2),"Padding","tight","TileSpacing","none");


% Plot Bars
ax1 = nexttile(uptl);cla;
b = bar(plot_means', 'grouped', 'EdgeColor', 'k','BarWidth', 1, 'FaceColor','flat','CData',cmap);
hold on;
ylabel(sprintf(upperYLabel));

% Set Bar Colors
for i = 1:length(b)
    b(i).FaceColor = "flat";
    b(i).CData = cmap(i,:);
end

% Add errorbars
for i = 1:num_conc
    x_coords = b(i).XEndPoints;
    errorbar(x_coords, plot_means(i,:), repmat(nan,1,length(plot_means(i,:))), plot_sems(i,:), 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize',2);
end
grid off;
set(ax1, 'XTickLabel', {});

% Make Annotations using Dunnett's Test
ax2 = nexttile(lwrtl);cla;
hold on;
axis off;

% Get Bar x-locations
bar_x_matrix = zeros(num_gt, num_conc);
for i = 1:num_conc
    bar_x_matrix(:, i) = b(i).XEndPoints';
end

% Annotation y-locations
y_genotype   = 0.5;
y_conc_label = 1.8;
y_dunnett    = 2.8;

offset = 5;
pw_dunnett_tracks = offset:1.1:(offset+1.1*(length(concentrations)-1));
pw_dunnett_colors = num2cell(cmap_annotations,2)';


% Genotype comparison labels
for c = 1:num_conc
    track_y = pw_dunnett_tracks(c);
    track_col = pw_dunnett_colors{c};

    mc = returnDunnettMultCompareStats(squeeze(data(:,c,:)));
    for i = 1:size(mc,1)
        if contains(p_to_stars(mc(i,end)),"*")
            ns_offset = -0.1;
        else 
            ns_offset = 0.1;
        end
        if i==2 
            offset=-0.2; 
        else 
            offset=0; 
        end
        if i==2 
            txtoffset=-0.6+ns_offset; 
        else 
            txtoffset=0.3+ns_offset; 
        end
        text(mean([bar_x_matrix(1,c),bar_x_matrix(i+1,c)]), track_y + txtoffset, p_to_stars(mc(i,end)), 'HorizontalAlignment', 'center', ...
        'Color', track_col,'FontWeight', 'normal');
        plot([bar_x_matrix(1,c),bar_x_matrix(i+1,c)], [track_y+offset, track_y+offset], '-', 'Color', track_col, 'LineWidth', 1);
    end
end

for g = 1:num_gt
    % Plot Dunnett's Test Results
    mc = returnDunnettMultCompareStats(squeeze(data(:,:,g)));
    stat_string = ["Ctrl"];
    for i = 1:size(mc,1)
        stat_string = [stat_string,p_to_stars(mc(i,end))];
    end
    text(bar_x_matrix(g,:)', repmat(y_conc_label,1,length(stat_string)), string(concentrations), 'HorizontalAlignment', 'center', 'FontSize', 9);
    text(bar_x_matrix(g,:)', repmat(y_dunnett,1,length(stat_string)), stat_string, 'HorizontalAlignment', 'center', 'FontWeight', 'normal', 'Color', 'b');

    % Genotype Spreader Bars
    gt_bars_x = bar_x_matrix(g, :);
    pad = 0.12;
    plot([min(gt_bars_x)-pad, max(gt_bars_x)+pad], [y_genotype+0.4, y_genotype+0.4], '-k', 'LineWidth', 1);
    text(mean(gt_bars_x), y_genotype, genotypes{g}, 'HorizontalAlignment', 'center', 'FontWeight', 'normal');
end


xlim(xlim(ax1)); % Match X-limits to the bar plot
ylim([0, max(pw_dunnett_tracks)+1.3]);  % Scale annotations area of plot

end

function stars = p_to_stars(p)
if p <= 0.001,     stars = '***';
elseif p <= 0.01,  stars = '**';
elseif p <= 0.05,  stars = '*';
else,              stars = 'ns';
end
end

function mc = returnDunnettMultCompareStats(data)
[p,anovatab,stats] = anova1(data(:,:,1),[],"off");
mc = multcompare(stats,"ControlGroup",1,"Approximate",1,"CriticalValueType","dunnett","Display","off");
end
