% Function to calculate return levels and produce a return level plot for a
% series dataset.

function [rl, datasort] = returnlevel(data,Y);

a = 0.44;
N = length(data);
datasort = sort(data);
datasort = datasort(end:-1:1);
rl = (Y - 2*a) ./ ((1:1:N)' - a);

plot(log(rl), datasort, '.-')
xlabel('n(years)')
ylabel('z(n)')
Title('Return Level Plot')
set(gca,'XTick',log([1;2;5;10;20;50;100;200;500]))
set(gca,'XTickLabel',{num2str([1;2;5;10;20;50;100;200;500])})
set(gca,'Xlim',[log(1) log(500)])
grid on