% function Hout = hist2d(Dx,Dy,Xbin,Ybin)
% sum(Hout(:)) will always equal length(Dx) as value below the minimum bin
% values are grouped into the minimum bin and those above the highest bin
% value grouped into the highest bin.
% 
% Calculates and returns the 2 Dimensional Histogram of D. 

% 
function H = hist2d_JN(Dx,Dy,Xbin,Ybin)

Xn = length(Xbin) ;
Yn = length(Ybin) ;

N = length(Dx) ;
if length(Dy)~=N, error('length(Dx)~=length(Dy)'), end
   
H = zeros(Yn,Xn) ;

for i = 1:N
    x = find(Xbin <= Dx(i),1,'last');
    y = find(Ybin <= Dy(i),1,'last');
    if isempty(x), x=1; end
    if isempty(y), y=1; end
        H(y,x) = H(y,x)+1;
end

if sum(H(:)) ~= length(Dx), error('not all records included in histogram'), end