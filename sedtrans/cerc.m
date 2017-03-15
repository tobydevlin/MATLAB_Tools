% Qsy = cerc(K,hb,Hb,Alpb,T,rhof,rhos,g,p,units)

function Qsy = cerc(K,hb,Hb,Alpb,T,rhof,rhos,g,p,units)

sgn = sign(Alpb);
Alpb = abs(Alpb);
switch units
    case 'm3/day'
        K2 = 86400 * K / ((rhos-rhof)*g*(1-p));
    case 'm3/annum'
        K2 = 365.25 * 86400 * K / ((rhos-rhof)*g*(1-p));
    otherwise
        error('units must be either ''m3/day'' or ''m3/annum''')
end
khb = kh_calc(4*pi^2/g./T.^2 .* hb);
Cgb = 0.5 * (g*T/(2*pi).*tanh(khb)) .* (1+2*khb./sinh(2*khb));
ECgb = 1/8*rhof*g*Hb.^2.*Cgb;
Qsy = sgn .* K2 .* ECgb .* sin(Alpb/180*pi) .* cos(Alpb/180*pi);