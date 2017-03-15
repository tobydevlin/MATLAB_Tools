% Function to calculate the inverse (or quantile) beta cumulative
% probability distribution function

function q = qbeta(p, alp, bet)

q = fminbnd(@fmin, 0, 1, [], alp, bet, p);

function err = fmin(q, alp, bet, p)
err = (pbeta(q, alp, bet) - p).^2;