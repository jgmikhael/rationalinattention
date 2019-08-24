% This code simulates the effect of dopamine on interval timing, both in 
% reproduction (Malapani et al., 1998; Lake & Meck, 2013) and estimation
% (Soares et al., 2016; Maricq et al., 1981), as well as the effect of
% controllability on clock speed (Bizo & White, 1995).
% Written 12Aug19 by JGM.

clear; close all; clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Malapani et al. (1998)

mu = [8 21];                    % timed durations
DA = [.2; 1]*[1 1];             % DA levels during [encoding, decoding]
t = 0:.1:35;                    % time domain
k0 = .05;                       % unit cost of information per time

[like, post] = deal(nan(length(t),4));
d = [2 4; 1 3]; % ordering for easy assignment to colors, legend
for e = 1:2
    [~,~ ,prior(:,e),like(:,d(e,:)),post(:,d(e,:))] = ...
        TimeModel(mu,DA(e,:),k0,t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(101)
for e = 1:2
    subplot(2,1,e)
    plot(t,prior(:,e),'k--')
    hold on
    h1 = plot(t,like(:,d(e,:)),'k');
    hold on
    h2 = plot(t,post(:,d(e,:)),'r');
    title(['DA = ' num2str(DA(e,2))])
    ylabel('Frequency')
    ylim([0 1.1*max(max(post))])
end
xlabel('Time (s)')
subplot(2,1,1)
legend('Prior','Likelihood, Short','Likelihood, Long',...
    'Posterior, Short', 'Posterior, Long')

figure(1)
figName{1} = 'MG98';
C = [.85 .25 .25; .3 .65 .65; .6 .2 .2; .3 .45 .45]; % color scheme
h = plot(t,post);
set(h, {'color'}, num2cell(C,2));
xlim([0 30])
xlabel('Time (s)')
ylabel('Frequency')
legend('8s High', '8s Low', '21s High', '21s Low')  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lake & Meck (2013)

muList = [7 17];             	% timed durations
[l, b, h] = deal(.8, 1, 1.2);  	% low, baseline, high DA
DA = [b h; b l; b b];         	% DA levels
t = 0:.1:40;                  	% time domain
k0 = .1;                        % unit cost of information per time

[like, post] = deal(nan(length(t),2,3));
for k = 1:2                    	% short, long duration
    mu = muList(k);
    for e = 1:3                	% low, baseline, high DA
        [~,~,~,~,post(:,k,e)] = TimeModel(mu, DA(e,:),k0,t);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
figName{2} = 'LM13';
C = [1 0 0; 0 0 1; 0 0 0];     	% color scheme
lines = {'--','-'};           	% solid vs dashed lines
for k = 1:2
    for e = 1:3
        p = post(:,k,e)./max(post(:,k,e));
        plot(t,p,'Color',C(e,:),'LineStyle',lines{k})
        hold on
    end
end
xlim([0 40])
xticks(5*(0:8))
yticks((0:.2:1))
xlabel('Time (s)')
ylabel('Proportion Max Response Rate')
legend('7s High','7s Low','7s Baseline','17s High','17s Low','17s Baseline')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bizo & White (1995)

Rt = .01:.01:.05;             	% total reinforcer rate
Rc = .01;                      	% rewards contingent on timing
e = .01;                        % smoothing term (epsilon)
R = log(Rt+e)-log(Rt-Rc+e)+e;   % added value of timing-contingent rewards
DA = R'*[1 1];                  % DA levels

k0 = .1;                        % unit cost of information per time
dur = 1;                        % timed duration (arbitrary)
t = 0:.01:dur;                  % time domain (arbitrary)

% compute scaling factor (i.e., pacemaker rate)
l = length(DA);
eta = nan(l,2);
for e = 1:l
   [eta(e,:),~,~,~,~] = TimeModel(dur, DA(e,:), k0, t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
figName{3} = 'BW95';
subplot(2,1,1)
plot(Rt,1./eta,'-ko')
xlabel('Obtained Reinforcer Rate (a.u.)')
xticks(.01:.01:.05)
ylabel('Pacemaker Period (a.u.)')

subplot(2,1,2)
plot(Rc./Rt,1./eta,'-ko')
xlabel('Proportion of Reinforcers Obtained From Timing')
xticks(.2:.2:1)
ylabel('Pacemaker Period (a.u.)')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soares et al. (2016)

dur = [.6, 2.4];              	% short and long durations
DA = [.7; .3]*[1 1];         	% DA levels
t = 0:.01:5;                    % time domain
k0 = .2;                        % unit cost of information per time

% compute p(Long)
p = nan(2,length(t));
for e = 1:2
    [~,~,~,~,post] = TimeModel(dur, DA(e,:), k0, t);
    p(e,:) = post(:,2)./sum(post,2);

	figure(401)
	subplot(1,2,e)
	plot(t,post)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
figName{4} = 'SP16';
plot(t,p(2,:),'r')
hold on
plot(t,p(1,:),'k')
xlim([dur(1)-.05,dur(2)+.05])
xticks([0.6 1.5 2.4])
yticks([0 0.5 1])
xlabel('Time Interval')
ylabel('p(Long Choice)')
legend('Low', 'High','Location','Northwest')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maricq et al. (1981)

pairs = [1 4; 2 8; 4 16];       % pairs of timed durations
[b, h] = deal(.7, .8);        	% baseline, high DA
DA = [b b; b h];             	% DA levels
t = 0:.1:20;                    % time domain
k0 = .2;                        % unit cost of information per time

% compute p(Long)
p = nan(3,length(t),2);
for e = 1:3                   	% for each interval pair
	for k = 1:2               	% for each DA pair
        [~,~,~,~,post] = TimeModel(pairs(e,:), DA(k,:), k0, t);
        p(e,:,k) = post(:,2)./sum(post,2);
        
      	figure(501)
        subplot(1,3,e)
        plot(t,post)
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
figName{5} = 'MC81';
C = [0 .3 .6]'*[1 1 1];         % color scheme
for e = 1:3
	int = t>pairs(e,1) & t<pairs(e,2);
	plot(t(int),p(e,int,1),'Color',C(e,:))
    hold on
    plot(t(int),p(e,int,2),'--','Color',C(e,:))
    hold on
end
xlabel('Signal Duration (s)')
ylabel('Proportion ''Long'' Response')
xlim([0 16])
ylim([0 1.1])
yticks([0 1])
legend('Baseline, 1s vs 4s','High, 1s vs 4s', 'Baseline, 2s vs 8s', ...
    'High, 2s vs 8s', 'Baseline, 4s vs 16s', 'High, 4s vs 16s',...
    'Location','Southeast')