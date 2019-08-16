function [  ] = plot_speed_through_generations( listGenes, number,generation)
importModelCodons();
%Plotting the change in speed through 100 generations
list_codon_times = zeros(length(listGenes),length(listGenes(1).Codons));
list_ratios = zeros(1,length(listGenes));
size(list_codon_times)
for i=1:length(listGenes)
    list_codon_times(i,:) = plot_speed(listGenes(i),number);
    list_ratios(i) = mean(list_codon_times(i,1:ceil(length(list_codon_times(i))/4))) / mean(list_codon_times(i,ceil(length(list_codon_times(i))/4):3*ceil(length(list_codon_times(i))/4)));
end

list_mean_times = mean(list_codon_times(:,:));

mean_beginning_time = mean(list_mean_times(1:(length(list_mean_times)/4)));
mean_end_time = mean(list_mean_times(3*(length(list_mean_times)/4):end));
mean_middle_time = mean(list_mean_times(2*(length(list_mean_times)/4):3*(length(list_mean_times)/4)));
mean_not_beginning_time = mean(list_mean_times((length(list_mean_times)/4):end))
beginning_end_ratio = mean_beginning_time/mean_end_time;
beginning_middle_ratio = mean_beginning_time/mean_end_time;
beginning_rest_ratio = mean_beginning_time/mean_not_beginning_time;

fileID = fopen('ratios.txt','a');
contents = strcat('At ',num2str(generation),', ', 'Gene ',num2str(number), ': Beginning: ', num2str(mean_beginning_time), ' Middle: ', num2str(mean_middle_time), ' End: ', num2str(mean_end_time), ' Beginning/End Ratio: ', num2str(beginning_rest_ratio));
fprintf(fileID,contents);
fclose(fileID);

plot(list_ratios)
title('Ratio of beginning vs end translation time vs Generation')
%plotting_function(list_mean_times,strcat('Average translation time vs position at generation ',num2str(generation),'for previous 100 generations. Total time: ',num2str(sum(list_mean_times))),number)
end