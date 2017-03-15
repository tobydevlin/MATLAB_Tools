% Qsy = kamphuis(K,Hb,Alpb,T,slope,d,rhos,p,units)

function Qsy = kamphuis(Hb,Alpb,T,slope,d,rhos,p,units)

sgn = sign(Alpb);
Alpb = abs(Alpb);
switch units
    case 'm3/day'
        K2 = 86400 * 2.27 / rhos / (1-p);
    case 'm3/annum'
        K2 = 365.25 * 86400 * 2.27 / rhos / (1-p);
    otherwise
        error('units must be either ''m3/day'' or ''m3/annum''')
end
Qsy = sgn .* K2 .* Hb.^2 .* T.^(1.5) .* slope.^(0.75)...
    .* (d).^(-0.25) .* (sin(2*Alpb/180*pi)).^0.6;