% This code illustrates the effects of information cost and prior precision
% on the relationship between reward incentives and likelihood precision.
% Written 2Jul20 by JGM.

clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
l0L = [0 .5 1 1.5 2];             	% list of lambda_0's (prior precision)
kL = [.1 .2 .3 .5 1];             	% list of kappas (unit cost of information)
R = 0:.001:.5;                      % reward incentive

figure(1); figName{1} = 'RvsLambda1';

% vary lambda_0
subplot(1,2,1)
C = linspecer(length(l0L));         % color scheme (available at MathWorks)
k = .5;                             % fixed kappa (unit cost of information) 
for e = length(l0L):-1:1
    l0 = l0L(e);                    % lambda_0 (prior precision)
    l = max(0, 2*R./k - l0);        % lambda (likelihood precision) 
    h(length(kL)-e+1) = plot(R, l,'Color',C(e,:));
    hold on
end
legend(fliplr(h),'\lambda_0 = 0', '\lambda_0 = 0.5', '\lambda_0 = 1', ...
    '\lambda_0 = 1.5', '\lambda_0 = 2','Location','Northwest',...
    'Box','Off','Interpreter', 'tex');

% vary kappa
subplot(1,2,2)
C = linspecer(length(kL));          % color scheme
l0 = .5;                            % fixed lambda_0
for e = length(kL):-1:1
    k = kL(e);                      % kappa
    l = max(0, 2*R./k - l0);        % lambda
    h(length(kL)-e+1) = plot(R, l,'Color',C(e,:));
    hold on
end
legend(fliplr(h),'\kappa = 0.1', '\kappa = 0.2', '\kappa = 0.3', ...
    '\kappa = 0.5','\kappa = 1','Location','Northwest',...
    'Box','Off','Interpreter', 'tex');

% ttl = {'\lambda_0','\kappa'};
for e = 1:2
    subplot(1,2,e)
    % title(['Changes in ' ttl{e}],'Interpreter', 'tex')
    xlabel('Reward Incentive\it R','Interpreter', 'tex')
    ylabel('Likelihood Precision \lambda*','Interpreter', 'tex')
    set(gca,'box','off')
end