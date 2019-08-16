function [  ] = plot_frequency( actual, random, name, number)
xvec = linspace(1,length(actual),length(actual)); % x vector for plotting
%plot(1:length(struct.a(:,1)),struct.a(:,1),1:length(struct.a(:,1)),struct.a(:,2));% Plots the actual and random distribution
plot(xvec,actual,xvec,random);
legend('Actual distribution', 'Random distribution');
set(gca,'xtick',1:(ceil(length(xvec)/30)):length(xvec),'xticklabel',1:(ceil(length(xvec)/30)):length(xvec));
xlabel('Position on Sequence')
ylabel('Frequency of Most Expressed Codon')
title('Frequency of Most Expressed Codon vs Position on Sequence')
savename = strcat('./Results/Plots/Frequency/gene',num2str(number),'/freq_',name,datestr(datetime('now')),'.fig');
savefig(gcf,savename)
end