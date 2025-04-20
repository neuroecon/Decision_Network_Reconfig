% prepare data for figures
%% Figure 1 - no data
%% Figure 2
% 2A: FC matrix
% Define cluster sizes and total matrix size
cluster_sizes = [35, 35, 20, 10];
n = sum(cluster_sizes);  % should be 100

% Pre-allocate matrix
M = zeros(n, n);

% Determine indices for each cluster
idx = cell(4,1);
start_idx = 1;
for k = 1:4
    idx{k} = start_idx:(start_idx + cluster_sizes(k) - 1);
    start_idx = start_idx + cluster_sizes(k);
end

% Fill matrix with baseline values according to specified rules
for i = 1:4
    for j = i:4  % Only fill upper triangular part
        if i == j
            % Within cluster: values in [0.6, 1]
            block = 0.6 + (1 - 0.6) * rand(length(idx{i}), length(idx{i}));
        else
            % Between clusters:
            if ( (i == 1 && j == 3) || (i == 1 && j == 4) || (i == 3 && j == 4) )
                % For these specific pairs: values in [0.1, 0.4]
                block = 0.1 + (0.4 - 0.1) * rand(length(idx{i}), length(idx{j}));
            else
                % For all other between-cluster pairs: values in [-0.4, -0.1]
                block = -0.4 + ((-0.1) - (-0.4)) * rand(length(idx{i}), length(idx{j}));
            end
        end
        % Place the block in the upper triangular section
        M(idx{i}, idx{j}) = block;
        % Ensure symmetry by mirroring the block to the lower triangular section
        M(idx{j}, idx{i}) = block';
    end
end


% Visualize the matrix using a heatmap
figure;
imagesc(M);
axis square;                     % Force the heatmap to be square
set(gca, 'XTick', [], 'YTick', []);  % Remove tick marks
set(gca, 'box', 'off');          % Remove the box around the axes

% Define hex colors and convert them to normalized RGB values
lower_color = [hex2dec('A0'), hex2dec('E7'), hex2dec('E5')] / 255; % for -1
upper_color = [hex2dec('FF'), hex2dec('AE'), hex2dec('BC')] / 255; % for 1
mid_color   = [1, 1, 1];  % white for 0

% Create a diverging colormap
nColors = 256;  % Total number of colors in the colormap
nHalf   = floor(nColors/2);

% Interpolate from lower_color to mid_color
cmap_lower = [linspace(lower_color(1), mid_color(1), nHalf)', ...
              linspace(lower_color(2), mid_color(2), nHalf)', ...
              linspace(lower_color(3), mid_color(3), nHalf)'];
          
% Interpolate from mid_color to upper_color
cmap_upper = [linspace(mid_color(1), upper_color(1), nColors - nHalf)', ...
              linspace(mid_color(2), upper_color(2), nColors - nHalf)', ...
              linspace(mid_color(3), upper_color(3), nColors - nHalf)'];

custom_cmap = [cmap_lower; cmap_upper];

% Apply the custom colormap and adjust color limits
colormap(custom_cmap);
caxis([-1 1]);  % Ensures that -1 maps to the lower_color and 1 to the upper_color

%% Figure 3 - in R

