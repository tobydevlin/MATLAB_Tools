% Function to create a custom figure with n axes downwards and the length
% and width specified in centimeters
% This makes report quality figures quickly
%
%   Inputs:
%           heigth - the height of the figure in cm
%           width  - the width of the figure in cm
%           n_axes - the number of axes downwards
% 
%   Output: 
%           f  - handle to the figure
%
%
% TDevlin Mar2016

function [ f ] = figureTD( width , height )

    f = figure('units','centimeters','paperunits','centimeters','position',[10 10 width height],'PaperPosition',[0 0 width height]);

end