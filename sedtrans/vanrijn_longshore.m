% Qsy = vanrijn_longshore(Hb,Alpb,T,slope,d,...
%                                    Vtide,rhos,g,p,units)
% From ICCE 2006

function Qsy = vanrijn_longshore(Hb,Alpb,T,slope,d,...
                                    Vtide,rhos,g,p,units)
sgn = sign(Alpb);
Alpb = abs(Alpb);

Kswell = ones(size(Hb));
Kswell(T>6) = T(T>6)/6;
Kgrain = (0.2/1000) ./ d;
Kgrain(Kgrain<0.1) = 0.1;
Kslope = slope / 0.01;
Kslope(Kslope>1.25) = 1.25;
Kslope(Kslope<0.75) = 0.75;
    
switch units
    case 'm3/day'
        K2 = 86400 * Kswell .* Kgrain .* Kslope ./ rhos ./ (1-p);
    case 'm3/annum'
        K2 = 365.25 * 86400 * Kswell .* Kgrain .* Kslope ./ rhos ./ (1-p);
    otherwise
        error('units must be either ''m3/day'' or ''m3/annum''')
end
Vwave = sgn .* 0.3 .* (g * Hb).^0.5 .* sin(2*Alpb);
Veff = Vwave + Vtide;
Qsy = 42 * K2 .* Hb.^2.5 .* Veff;