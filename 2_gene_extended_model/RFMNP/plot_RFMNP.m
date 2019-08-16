function plot_RFMNP(T, Rate, Density, Pool, Init, mRNA)
% plot_RFMNP(T, Rate, Density, Pool)
%   animate the ribosome densities on multiple genes through time.
%
% June 2015
% Alon Diament / Tuller Lab

nORFs = length(Density);
xlims = [1, max(cellfun(@(x) size(x, 1), Density))];
color_map = colormap;

Translating = sum(cell2mat(cellfun(@(x, y) y*sum(x(1:end-1, :), 1), Density, ...
    mat2cell(mRNA, ones(nORFs, 1)), 'UniformOutput', false)), 1); % ribosomes

for t = 1:length(T)
    min_init = Inf;
    max_init = -Inf;
    for o = 1:nORFs
        G = Init(o) .* tanh(Pool(1:t));
        this_color = round(o / nORFs * size(color_map, 1));
        this_color = color_map(this_color, :);
        subplot(2, 1, 1);
        plot(Density{o}(1:end-1, t), 'o-', 'MarkerFaceColor', this_color, ...
            'Color', this_color);
        hold on
        if o == nORFs
            hold off
        end
        xlim(xlims);
        xlabel('chunk');
        ylabel('ribosomal density (x_i)');
        subplot(2, 3, 4);
        plot(T(1:t), Rate{o}(1:t), '.-', 'MarkerFaceColor', this_color, ...
            'Color', this_color);
        hold on
        if o == nORFs
            hold off
        end
        xlabel('time');
        ylabel('translation rate');
        subplot(2, 3, 5);
        init_vec = G.*(1-Density{o}(1, 1:t));
        plot(T(1:t), init_vec, '.-', 'MarkerFaceColor', this_color, ...
            'Color', this_color);
        tframe = max(1, t-30):t;
        % zoom in on time-frame
        min_init = min(min_init, min(init_vec(tframe)));
        max_init = max(max_init, max(init_vec(tframe)));
        hold on
        if o == nORFs
            hold off
        end
        xlabel('time');
        ylabel('initiation rate');
    end
    subplot(2, 3, 6);
    plot(T(1:t), Pool(1:t), '.-', 'MarkerFaceColor', 0.5*[0,1,0], ...
        'Color', 0.5*[0,1,0]);
    hold on
    plot(T(1:t), Translating(1:t), '.-', 'MarkerFaceColor', 0.5*[1,0,0], ...
        'Color', 0.5*[1,0,0]);
    hold off
    xlabel('time');
    ylabel('ribosomes');
    subplot(2, 3, 5);
    ylim([min_init - eps, max_init + eps]);
    pause(0.1);
end
% legends
subplot(2, 3, 6);
legend('in pool', 'translating', 'Location', 'Best');
subplot(2, 1, 1);
legend(arrayfun(@(x) sprintf('gene-%02d', x), [1:nORFs]', ...
    'UniformOutput', false), ...
    'Location', 'NorthEastOutSide');
