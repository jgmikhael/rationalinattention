% This code simulates the effect of dopamine on reinforcement learning,
% for both the modulation of learning from positive and negative feedback
% (Cools et al., 2009) and control of the exploration-exploitation balance
% (Cinotti et al., 2019).
% Written 7Aug19 by JGM.

clear; close all; clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cools et al. (2009)

% parameters
DA = .4:.02:.6;                 % DA levels
k = .1;                         % unit cost of information
[mu0, R] = deal(DA);          	% prior mean, average reward
l0 = 5;                         % prior precision (arbitrary)
l = max(0,2*R/k-l0);           	% likelihood precision
mu = [0; 1];                    % likelihood means
r = mu(1)-1:.01:mu(2)+1;        % reward domain

% compute posteriors (after just 1 learning step)
lh = l+l0;                      % posterior likelihoods
muL = mu*ones(1,length(DA));
mu0L = [1; 1]*mu0;
muh = (l.*muL+l0*mu0L)./lh;     % posterior means

% distributions
prior = nan(length(r),length(DA));
for k = 1:length(DA)
    prior(:,k) = normpdf(r,mu0(k),1/sqrt(l0));
end

like = nan(length(r),2,length(DA));
for k = 1:length(DA)
    for e = 1:2
        like(:,e,k) = normpdf(r,mu(e),1/sqrt(l(k)));
    end
end

[post, postC] = deal(nan(length(r),2,length(DA)));
for k = 1:length(DA)
    for e = 1:2
        post(:,e,k) = normpdf(r,muh(e,k),1/sqrt(lh(k)));
        postC(:,e,k) = normcdf(r,muh(e,k),1/sqrt(lh(k))); % cumulative
    end
end

% relative accuracy
ind = find(r==0.5);             % index for which r = 0.5
a = nan(length(DA),1);
for k = 1:length(DA)
    a(k) = (1-postC(ind,2,k)) - postC(ind,1,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(101)
for k = 1:length(DA)
    subplot(length(DA),1,k)
    plot(r,prior(:,k),'Color',.8*[1 1 1])
    hold on
    plot(r,like(:,:,k),'k--')
    hold on
    plot(r,post(:,:,k),'b')
end
subplot(length(DA),1,length(DA))
xlabel('Rewards')

figure(1)
figName{1} = 'CD09';
plot(DA,0*DA,'Color',.5*[1 1 1],'LineWidth',2)
hold on
scatter(DA,a,'k','filled')
xlim([min(DA) max(DA)])
xlabel('DA Level')
ylabel('Relative Learning')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cinotti et al. (2019)

L = [87.5 6.25 6.25]/100;      	% reward prob in low-risk condition
H = [62.5 18.75 18.75]/100;   	% reward prob in high-risk condition
arm = [L; H];

% parameters
antag = 0:.1:.3;             	% DA antagonist level
DAb = .9;                     	% baseline DA (arbitrary)
DA = DAb-antag';                % DA levels
R = DA;                        	% average reward
k = .1;                       	% unit cost of information
r = -1:.01:2;                  	% reward domain
beta = 4;                      	% inverse temperature for softmax

% prior
m0 = mean(L);                	% mean
l0 = 10;                     	% precision (arbitrary)

cond = {'Low-Risk','High-Risk'};% subfigure labels
for q = 1:2
    
    % likelihood
    m = arm(q,:);             	% mean
    l = max(0,2*R/k-l0);       	% precision
    
    % posterior
    mh = (l*m+l0.*m0)./(l+l0);	% mean
    lh = l+l0;                 	% precision
    
    % p(selecting one of two smaller rewards), according to softmax
    p(:,q) = 1./(1+2*exp(beta*(mh(:,1)-mh(:,2))));
    
    % distributions: prior, likelihood large/small, posterior large/small
    prior = normpdf(r,m0,1./sqrt(l0));
    [likeL,postL,likeS,postS] = deal(nan(length(r),2));
    for e = 1:length(DA)
        likeL(:,e) = normpdf(r,m(1),1./sqrt(l(e)))';
        postL(:,e) = normpdf(r,mh(e,1),1./sqrt(lh(e)))';
        likeS(:,e) = normpdf(r,m(2),1./sqrt(l(e)))';
        postS(:,e) = normpdf(r,mh(e,2),1./sqrt(lh(e)))';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(201)
    subplot(2,1,q)
    
    plot(r,prior,'k--')
    hold on
    h1 = plot(r,likeS);
    hold on
    h2 = plot(r,likeL);
    hold on
    h3 = plot(r,postS);
    hold on
    h4 = plot(r,postL);
    xlabel('Rewards')
    ylabel('Distribution')
    title(cond{q})
    
    v = linspace(.2,.8,length(DA))';
    C1 = v*[1 1 1];
    C2 = v*[.1 .1 .1]*5;
    C3 = v*[.1 1 .1];
    C4 = v*[.1 .1 1];
    
    set(h1, {'color'}, num2cell(C1,2));
    set(h2, {'color'}, num2cell(C2,2));
    set(h3, {'color'}, num2cell(C3,2));
    set(h4, {'color'}, num2cell(C4,2)); 
end

figure(2)
figName{2} = 'CK19';
b = bar(antag,100*p,'k');

c = [62 84 141; 143 51 52]/255; % color scheme
for q = 1:2
    b(q).FaceColor = 'flat';
    b(q).EdgeColor = 'flat';
    b(q).CData = c(q,:);
end

ylim([0 40])
xlabel('DA Antagonist')
ylabel('% Exploration')
xticks(antag)
yticks(0:10:40)