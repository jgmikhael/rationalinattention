% This code illustrates the effect of likelihood precision on the posterior
% distribution for the case of two reward sources.
% Written 1Aug19 by JGM.

clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
figName{1} = 'Bayes';

s = 4; l = 8;                   % small and large reward
r = 0:.1:12;                    % all possible reward values
sigL = [.8, 2];                 % low and high SD (high and low precisions)

p = nan(2,length(r));
maxY = 0;
for e = 1:2
    sig = sigL(e);
    
    % set likelihoods p(s|t), p(l|t)
    ps = normpdf(r,s,sig); ps = ps./sum(ps);
    pl = normpdf(r,l,sig); pl = pl./sum(pl);
    
    % compute prior p(t)
    params = fitdist([s; l],'Normal');
    q = params.mu;
    x = params.sigma;
    pr = normpdf(r,q,x); pr = pr./sum(pr);
    
    % compute posteriors
    pos = ps.*pr; pos = pos./sum(pos);
    pol = pl.*pr; pol = pol./sum(pol);
    
    % find posterior means
    [~,v1] = max(pos); [~,v2] = max(pol);
    mps = r(v1); mpl = r(v2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(1,2,e)
    
    bl = -.002;                     % for visualization of blue line
    C = [0 .6 0; .8 0 0];           % color scheme
    titles = {'High Precision','Low Precision'};
    plot(r,pr,'--k','LineWidth',3)
    hold on
    plot(r,ps,'--','Color', .7*[1 1 1],'LineWidth',3)
    hold on
    plot(r,pl,'--','Color', .5*[1 1 1],'LineWidth',3)
    hold on
    plot(r,pos,'-','Color',C(e,:),'LineWidth',4)
    hold on
    plot(r,pol,'-','Color',C(e,:)/1.5,'LineWidth',4)
    hold on
    plot([mps mpl], bl+[0 0],'b') 	% highlight difference in rewards
   
    xlim([0 max(r)])
    xlabel('Reward')
    ylabel('Probability')
    title(titles{e})
    maxY1 = max([pos pol]);
    maxY = max([maxY1 maxY]);
    ylim([bl 1.03*maxY])
end
legend('Prior','Likelihood, Small','Likelihood, Large',...
        'Posterior, Small','Posterior, Large')