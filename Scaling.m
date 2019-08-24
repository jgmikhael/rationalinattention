% This code illustrates the effect of dopamine on time cell response
% profiles, when plotted against objective time.
% Written 1Aug19 by JGM.

clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
n = 20;                 % number of cells
DA = [1 1.5];           % [low, high] DA
k = .05;                % unit cost of information per second
m = 1:n;                % means, in subjective time
a = 1;                  % SD, in subjective time (arbitrary)
eta = a*sqrt(2*DA/k);   % two example eta's (scaling factors)       
mu = exp(m'./eta);      % means, in objective time
sig = a*mu./eta;        % SDs, in objective time
t = 0:.01:max(max(mu)); % objective time

% time cells: x(i,j) = (objective time i, cell j)
xLow = normpdf(t', mu(:,1)', sig(:,1)');  xLow = xLow./max(xLow);
xHigh = normpdf(t', mu(:,2)', sig(:,2)'); xHigh = xHigh./max(xHigh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
figName{1} = 'scaling';

C = linspecer(n);       % color scheme (available at MathWorks)
sp = {'                                                 '};

subplot(2,1,1)
h = plot(t-1,xLow);
set(h, {'color'}, num2cell(C,2));
title(strcat({'Low DA'}, sp))

subplot(2,1,2)
h = plot(t-1,xHigh);
set(h, {'color'}, num2cell(C,2));
title(strcat({'High DA'}, sp))

for e = 1:2
    subplot(2,1,e)
    ylabel('Time Cell Activation')
    set(gca,'ytick',[])
    xlim([0 5])
    ylim([0, max(xLow(:))])
    xlabel('Objective Time')
end