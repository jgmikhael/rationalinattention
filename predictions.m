% This code simulates predictions of the rational inattention framework.
% Written 2Mar20 by JGM.

clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% high controllability vs low controllability
m = [5 10];                     % likelihood means for low/high reward
l = [1 .1];                     % likelihood precision for high/low control
m0 = mean(m);                   % prior mean
params = fitdist(m','Normal');
l0 = 1/(params.sigma)^2;        % prior precision (approximate)
lp = l+l0;                      % posterior precision
range = 0:.1:1.5*m(2);

% high, low controllability
for c = 1:2
    k = l(c)./(l(c)+l0);
    mA(c) = k*m(1) + (1-k)*m0;
    mB(c) = k*m(2) + (1-k)*m0;
end

A1 = normpdf(range, mA(1), 1/sqrt(lp(1)));
A2 = normpdf(range, mA(2), 1/sqrt(lp(2)));

B1 = normpdf(range, mB(1), 1/sqrt(lp(1)));
B2 = normpdf(range, mB(2), 1/sqrt(lp(2)));

% softmax
T = 1;                              % temperature parameter
pA1 = 1/(1+exp(-(mA(1)-mA(2))/T));  % A1 vs A2
pB1 = 1/(1+exp(-(mB(1)-mB(2))/T));  % B1 vs B2

pA2 = 1-pA1;
pB2 = 1-pB1;



% for Prediction 3

mm = mean(m); % medium reward
for c = 1:2
    k = l(c)./lp(c);
    mC1(c) = k*mm+(1-k)*mean([m(1) mm]);
    mC2(c) = k*mm+(1-k)*mean([mm m(2)]);
    mAx(c) = k*m(1)+(1-k)*mean([m(1) mm]);
    mBx(c) = k*m(2)+(1-k)*mean([mm m(2)]);
end
C1 = normpdf(range, mC1(1), 1/sqrt(lp(1)));
C2 = normpdf(range, mC2(1), 1/sqrt(lp(1)));


Ax1 = normpdf(range, mAx(1), 1/sqrt(lp(1)));
Ax2 = normpdf(range, mAx(2), 1/sqrt(lp(2)));

Bx1 = normpdf(range, mBx(1), 1/sqrt(lp(1)));
Bx2 = normpdf(range, mBx(2), 1/sqrt(lp(2)));


Tg = [.1 1]; % temperature parameter for high, low gain control

% (a) gain control
pAhi_gain = 1/(1+exp(-(mC2(1)-mC1(1))/Tg(1)));
pAlo_gain = 1/(1+exp(-(mC2(1)-mC1(1))/Tg(2)));

% (b) rational inattention
pAhi_RI = 1/(1+exp(-(mC2(1)-mC1(1))/T));
pAlo_RI = 1/(1+exp(-(mC2(2)-mC1(2))/T));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); figName{1} = 'predictions1';
col = get(gca,'ColorOrder');

subplot(2,3,1)
plot(range,A1)
hold on
plot(range,B1)
hold on
plot(range,A2,'--','Color',col(1,:))
hold on
plot(range,B2,'--','Color',col(2,:))
xlabel('Reward')
ylabel('Probability')
yticks([0 .5])
ylim([0 .8])
legend('High Controllability, Small (A1)', ...
    'High Controllability, Large (B1)', ...
    'Low Controllability, Small (A2)', ...
    'Low Controllability, Large (B2)','Location','Northeast')
legend('A1','B1','A2','B2','Location','Northeast')
legend boxoff

subplot(2,3,4)
bar(1,pA1-pA2)
hold on
bar(2,pB1-pB2)
xticks([1 2])
xticklabels({'A','B'})
ylabel('Probability Difference')
% yticks(0)
ylim([-.5 1])
hold on
plot([0 3],[0 0],'k','LineWidth',1)
xlim([0 3])
xlabel('Arm')
legend('p(A1)-p(A2)', 'p(B1)-p(B2)', 'Location', 'Northeast')
legend boxoff


subplot(2,3,2)
plot(range,A1,'Color',col(1,:))
hold on
plot(range,A2,'Color',col(2,:))
hold on
plot(range,B1,'Color',col(1,:))
hold on
plot(range,B2,'Color',col(2,:))
xlabel('Reward')
ylabel('Probability')
ylim([0 .8])
legend('High Controllability','Low Controllability','Location','Northeast')
legend boxoff

subplot(2,3,5)
plot([1 2],[mA(1) mB(1)]-m0)
hold on
plot([1 2],[mA(2) mB(2)]-m0)
hold on
scatter([1 2],[mA(1) mB(1)]-m0,[],col(1,:),'filled')
hold on
scatter([1 2],[mA(2) mB(2)]-m0,[],col(2,:),'filled')
hold on
xticks([1 2])
xticklabels({'Small','Large'})
xlim([.5 2.5])
xlabel('Received Reward')
ylabel('Phasic DA')
ylim([-2.5 4])
plot([0 3],[0 0],'k','LineWidth',1)
legend('High Controllability','Low Controllability','Location','Northeast')
legend boxoff

subplot(2,3,3)
plot(range,Ax1,'--','Color',.5*[1 1 1])
hold on
plot(range,C1,'Color',col(1,:))
hold on
plot(range,C2,'Color',col(2,:))
hold on
plot(range,Bx1,'--','Color',.8*[1 1 1])
xlabel('Reward')
ylabel('Probability')
ylim([0 .8])
legend('A1','C1','C2','B2','Location','Northeast')
legend boxoff


subplot(2,3,6)
plot([1 2],[pAlo_gain pAhi_gain])
hold on
plot([1 2],[pAlo_RI pAhi_RI])
hold on
scatter([1 2],[pAlo_gain pAhi_gain],[],col(1,:),'filled')
hold on
scatter([1 2],[pAlo_RI pAhi_RI],[],col(2,:),'filled')
hold on
xticks([1 2])
xticklabels({'Low','High'})
xlabel('DA Level')
ylabel('Probability of Arm C2')
plot([0 3],[.5 .5],'k','LineWidth',1)
yticks([0 .5 1 1.5])
ylim([0 1.5])
xlim([.5 2.5])
legend('Gain Control', 'Rational Inattention', 'Location', 'Northeast')
legend boxoff


for e = 1:3
    subplot(2,3,e)
    yticks([0 .5])
end