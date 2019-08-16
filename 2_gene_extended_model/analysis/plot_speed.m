function [] = plot_speed(list_speed,number,generation)

%Declare a basic date format
date_format = 'mm_dd_yy';
%Make our figure to be saved later
figure(2)
%Set the title and labels
title(strcat('Translational time for each gene on sequence number ',num2str(number), ' Total time: ', num2str(sum(list_speed))))
xlabel('Position on Sequence')
ylabel('Time')
hold on
%Plot the ratios, transposing first.
list_speed = list_speed.';
plot(list_speed)

hold off
%xticks(custom_xticks)
%xticklabels(string(custom_xticks))
savename = strcat('Plots/Speed/',num2str(number),'_speed_',num2str(generation),'_',datestr(datetime('now'),date_format),'.fig');
savefig(savename)
end
