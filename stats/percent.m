
function xp = percent(x,p)

if p<0 | p>1
    error('p must be between 0 and 1')
end

x = sort(x);
n = numel(x);
i = n*p;

