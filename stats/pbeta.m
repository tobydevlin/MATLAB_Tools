% Function to calculate the beta cumulative probability distribution
% function

function p = pbeta(q, alp, bet)

p = betainc(q, alp, bet);