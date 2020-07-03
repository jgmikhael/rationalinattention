% This code simulates the effect of dopamine on interval timing, both in 
% reproduction (Malapani et al., 1998; Lake & Meck, 2013) and estimation
% (Soares et al., 2016; Maricq et al., 1981), the effect of controllability
% on clock speed (Bizo & White, 1995) and post-reward delays (Blanchard
% et al., 2013), the effect of average reward on clock speed (Killeen &
% Fetterman, 1988), and the effect of prefeeding on the central tendency
% (Ward & Odum, 2007).
% Written 12Aug19 by JGM.
% Updated 17Mar20 by JGM, with post-reward delay simulation.
% Updated 29Jun20 by JGM, with average reward and prefeeding simulations.

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

Rt = 1:5;                       % total reinforcer rate
Rc = 1;                      	% rewards contingent on timing
R = Rc./Rt;                     % value of timing-contingent rewards
DA = R'*[1 1];                  % DA levels

k0 = .1;                        % unit cost of information per time
dur = 1;                        % timed duration (arbitrary)
t = 0:1:dur;                    % time domain (arbitrary)

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
xticks(1:5)

subplot(2,1,2)
plot(Rc./Rt,1./eta,'-ko')
xlabel('Proportion of Reinforcers Obtained From Timing')
xticks(.2:.2:1)

for e = 1:2
    subplot(2,1,e)   
    yticks(.2:.1:.6)
    ylim([.2 .6])
    ylabel('Pacemaker Period (a.u.)')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blanchard et al. (2013)

% compute scaling factor (i.e., pacemaker rate)
R = [1; .7];                    % pre, post-reward
DA = R*[1 1];                   % DA level
dur_pre = 1;                    % standard unit of pre-reward duration
dur = [0 1 2 3 4 5 6 10];       % timed durations for post-reward
t = 0:.01:dur;                  % time domain (arbitrary)

l = length(dur);
etaPRE = nan(l,2);
etaPOST = nan(l,2);

for e = 1:l
   [etaPRE(e,:),mhpre(e),~,~,~] = TimeModel(dur_pre, DA(1,:), k0, t);
   [etaPOST(e,:),mhpost(e),~,~,~] = TimeModel(dur(e), DA(2,:), k0, t);
end
mhStandard = mean(mhpost(1:7));

w = mhpost./mhpre;
wStand = mhStandard./mhpre(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
figName{6} = 'BH13';
bar(3:9,w([1 2 3 4 5 6 8]),'FaceColor',[204 205 207]/255)
hold on
bar(1,wStand,'FaceColor',[70 70 72]/255)
xticks([1 3:9])
xticklabels({'standard task','0','1','2','3','4','5','10'})
xlabel('buffer duration')
ylabel('w term')

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Killeen & Fetterman (1988)

R = 0:.5:10;                    % reinforcer rate
DA = R'*[1 1];                  % DA levels

k0 = .2;                        % unit cost of information per time
dur = 1;                        % timed duration (arbitrary)
t = 0:1:dur;                    % time domain (arbitrary)

% compute scaling factor (i.e., pacemaker rate)
l = length(DA);
eta = nan(l,2);
for e = 1:l
   [eta(e,:),~,~,~,~] = TimeModel(dur, DA(e,:), k0, t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7)
figName{7} = 'KF88';
plot(R,eta,'-ko')
xlabel('Reinforcement Density')
xticks(0:10)
ylabel('Clock Speed')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ward & Odum (2007)

dur = [2 8];                    % short and long durations
DA = [.7; .3]*[1 1];         	% DA levels
t = 0:.01:8;                    % time domain
k0 = .2;                        % unit cost of information per time

% compute p(Long)
p = nan(2,length(t));
for e = 1:2
    [~,~,~,~,post] = TimeModel(dur, DA(e,:), k0, t);
    p(e,:) = post(:,2)./sum(post,2);

	figure(801)
	subplot(1,2,e)
	plot(t,post)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8)
figName{8} = 'WO07';
plot(t,p(1,:),'k')
hold on
plot(t,p(2,:),'k--')
xlim([dur(1)-.05,dur(2)+.05])
xticks([2 5 8])
yticks(0:.2:1)
xlabel('Time Interval')
ylabel('p(Long Choice)')
legend('Baseline', 'Disrupt','Location','Northwest')
