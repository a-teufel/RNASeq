function [] = plotting_function(times_list,titlestr,number,generation,varargin)
stem(1:length(times_list),times_list);% Plots the translation times

legend('Codon Translation Times');
%set(gca,'xtick',1:(length(codonTimes)),'xticklabel',Genes(number).Codons);
ylabel('Translation Time')
title(titlestr)

disp(ceil(length(times_list)))
% Set the X-Tick locations so that every other gene is labeled.
Xt = 1:ceil(length(times_list)/3);
Xl = [1 ceil(length(times_list))];
disp(Xl)
set(gca,'XTick',Xt,'XLim',Xl);
% Add the codons as tick labels.
ax = axis;    % Current axis limits
axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4);  % Y-axis limits
% Place the text labels
if nargin == 5
t = text(Xt,Yl(1)*ones(1,length(Xt)),varargin{1}.Codons(1,1:ceil(length(varargin{1}.Codons)/3)));
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
end


savename = strcat('./Results/Plots/Speed/gene',num2str(number),'/speed_',num2str(generation),datestr(datetime('now')),'.fig');
savefig(savename)
end

