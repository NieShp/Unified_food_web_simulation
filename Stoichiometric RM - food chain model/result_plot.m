clear
load U_new
clf
subplot(122)
U = sortrows(U, 1, 'descend');
U = sortrows(U, 2);
U(:,1) = log10(U(:,1));
consumer = reshape(U(:,4),[length(unique(U(:,1))) length(unique(U(:,2)))]);
coexistence = reshape(U(:,7),[length(unique(U(:,1))) length(unique(U(:,2)))]);
coexistence(coexistence>1e-13)=1;
coexistence(coexistence<1e-13)=0;
coexistence(consumer<1e-13)=-1;
consumer(coexistence~=0) = 0;
consumer(coexistence~=0) = NaN;
imagesc(consumer);
colormap([1 1 1; parula(256)]);  
colorbar;
y_ticks = unique(U(:,1));y_ticks = sort(y_ticks, 'descend');
x_ticks = unique(U(:,2));
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', 1:30:length(y_ticks), 'YTickLabel', y_ticks(1:30:length(x_ticks))); 
set(gca, 'XTick', 1:20:length(x_ticks), 'XTickLabel', x_ticks(1:20:length(x_ticks)));
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 1);

subplot(121)
resource = reshape(U(:,3),[length(unique(U(:,1))) length(unique(U(:,2)))]);
resource(isnan(consumer))= 0;
imagesc(resource);
colormap([1 1 1; parula(256)]);  % 1 1 1 是白色，parula 是数据的颜色映射
colorbar;
y_ticks = unique(U(:,1));y_ticks = sort(y_ticks, 'descend');
x_ticks = unique(U(:,2));
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', 1:30:length(y_ticks), 'YTickLabel', y_ticks(1:30:length(x_ticks))); 
set(gca, 'XTick', 1:20:length(x_ticks), 'XTickLabel', x_ticks(1:20:length(x_ticks)));
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 1);
set(gcf, 'Position', [300, 300, 1600, 600]); 
% exportgraphics(gcf, 'Effect_of_Stoichiometry_on_biomass.png', 'Resolution', 300); 