function [] = plot_ratios(list_ratios,number,generation)

%Declare a basic date format
date_format = 'mm_dd_yy';
%Make our figure to be saved later
figure(1)
%Set the title and labels
title(strcat('Ratio of mean translation time of codon ramp to rest of sequence for gene ',num2str(number)))
xlabel('Generations')
ylabel('Ratios')
hold on
%Plot the ratios, transposing first.
list_ratios = list_ratios.';
plot(list_ratios)

hold off
legend('beginning/middle','beginning/end','beginning/all')
%xticks(custom_xticks)
%xticklabels(string(custom_xticks))
savename = strcat('Plots/Ratios/',num2str(number),'_ratio_',num2str(generation),'_',datestr(datetime('now'),date_format),'.fig');
savefig(savename)
end
