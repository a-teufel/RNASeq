function [] = plot_score(scores_list,number,generation)

%Declare a basic date format
date_format = 'mm_dd_yy';
%Make our figure to be saved later
figure(2)
%Set the title and labels
title(strcat('Network Translational Efficiency vs Generation, from 1 to ',num2str(generation)))
xlabel('Generation')
ylabel('Score')
hold on
%Plot the ratios, transposing first.
plot(scores_list)

hold off
%xticks(custom_xticks)
%xticklabels(string(custom_xticks))
savename = strcat('Plots/',num2str(number),'_score_',num2str(generation),'_',datestr(datetime('now'),date_format),'.fig');
savefig(savename)
end
