function [  ] = plot_speed( listGenes, name, number)
importModelCodons();
beginning_mean_times = zeros(1,length(listGenes));
for p=1:length(listGenes)
    codonTimes = zeros(1,length(listGenes(p).Codons)/4);
    for k =1:length(listGenes(p).Codons)/4 %for each codon in the list 
        Gene = listGenes(p);
        for i=1:length(AAA)
            if strcmp(Gene.Codons(k),AAA(i))
                codonTimes(k) = VarName2(i);
                %disp('in here');
            end
        end
    end
    beginning_mean_times(p) = mean(codonTimes);   
end


end_mean_times = zeros(1,length(listGenes));
for p=1:length(listGenes)
   end_codon_times = zeros(1,length(listGenes(p).Codons)/4);
    for k =(length(listGenes(p).Codons)/4):length(listGenes(p).Codons) %for each codon in the list 
        Gene = listGenes(p);
        for i=1:length(AAA)
            if strcmp(Gene.Codons(k),AAA(i))
                end_codon_times(k) = VarName2(i);
                %disp('in here');
            end
        end
    end
    end_mean_times(p) = mean(end_codon_times);   
end
total_mean_times = horzcat(beginning_mean_times, end_mean_times);
stem(1:length(codonTimes),codonTimes);% Plots the actual and random distribution
avg_beginningTime = sum(codonTimes(1:ceil(length(codonTimes)/4)))/length(codonTimes(1:ceil(length(codonTimes)/4)));
avg_endTime = sum(codonTimes(ceil(1*length(codonTimes)/4):end))/length(codonTimes(ceil(1*length(codonTimes)/4):end));
ratio = avg_beginningTime/avg_endTime;
legend('Codon Translation Times');
%set(gca,'xtick',1:(length(codonTimes)),'xticklabel',Genes(number).Codons);
ylabel('Translation Time')
plot_title = strcat('Time vs Position on Sequence ',num2str(number),' at 0. Total: ',num2str(sum(codonTimes)));
title(plot_title)


% Set the X-Tick locations so that every other gene is labeled.
Xt = 1:ceil(length(Genes.Codons)/2);
Xl = [1 ceil(length(Genes.Codons)/2)];
set(gca,'XTick',Xt,'XLim',Xl);
% Add the codons as tick labels.
ax = axis;    % Current axis limits
axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4);  % Y-axis limits
% Place the text labels
t = text(Xt,Yl(1)*ones(1,length(Xt)),Genes.Codons(1,1:ceil(length(Genes.Codons)/2)));
disp(t)
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'Rotation',90);
% Remove the default labels
set(gca,'XTickLabel','')
% Get the Extent of each text object.  This
% loop is unavoidable.
for i = 1:length(t)
  ext(i,:) = get(t(i),'Extent');
end
% Determine the lowest point.  The X-label will be
% placed so that the top is aligned with this point.
LowYPoint = min(ext(:,2));
% Place the axis label at this point
XMidPoint = Xl(1)+abs(diff(Xl))/2;
tl = text(XMidPoint,LowYPoint,'Position on Sequence', ...
          'VerticalAlignment','top', ...
          'HorizontalAlignment','center');
savename = strcat('./Results/Plots/Speed/gene',num2str(number),'/speed_',name,'_47k_',datestr(datetime('now')),'.fig');
fileID = fopen('avg_ratios.txt','a');
contents = strcat('At 0, Gene',num2str(number), ': Beginning: ', num2str(avg_beginningTime), ' Middle: ', num2str(avg_endTime), ' Ratio: ', num2str(ratio));
fprintf(fileID,contents);
fclose(fileID);
%savefig(savename)

end
