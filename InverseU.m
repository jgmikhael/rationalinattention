% This code illustrates the inverse-U-shaped relationship between dopamine
% and performance when a limit to true precision, but not estimated
% precision, is imposed.
% Written 14Aug19 by JGM.

clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Free parameters
l0 = 40;                                    % prior precision
L = 20;                                     % limit to true precision
d = 1:.01:100;                              % range of d

[~, dind] = min(abs(d-L));
before = 1:dind;
after = dind+1:length(d);

% error as a function of d
Ea = 1./(l0+d);                             % d < L
Eb = (l0*L+d.^2)./(L*(l0+d).^2);            % d > L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
figName{1} = 'inverseU';

plot(d(before), Ea(before),'k')
hold on
plot(d(after), Eb(after))
hold on
plot(d, Ea,'k--')
hold on
plot(d,Eb,'--','Color',[0 .447 .741])
hold on
plot(L,0,'ro')

set(groot,'defaultLegendInterpreter','LaTeX');
legend('$E(d)$ when $d \leq L$','$E(d)$ when $d > L$')
set(groot,'defaultLegendInterpreter','none');
xlabel('Desired Precision')
ylabel('Error')