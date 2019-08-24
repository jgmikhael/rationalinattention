% This function simulates the rational inattention model of interval
% timing.
% Written 12Aug19 by JGM.

function [eta, mh, prior, like, post] = TimeModel(mu, DA, k0, t)

% inputs
%   mu:         interleaved values                              (1xn)
%   DA:         DA at encode and decode                         (1x2)
%   k0:         unit cost of information per time               (1x1)
%   t:          time domain                                     (1xl)

% outputs
%   eta:        scaling factor                                  (1x1)
%   mh:         posterior means                                 (nx1)
%   prior:      prior distribution                              (lx1)
%   like:       likelihood distributions                        (lxn)
%   post:     	posterior distributions                         (lx1)

% parameters
R0 = DA;                            % average reward
a = 1;                              % SD in subjective space (arbitrary)

% retrieve eta from DA
eta = a*sqrt(2*R0/k0);           	% approximation from l = max(0,2R/k-l0)

% likelihoods (subjective)
m = log(mu+1)'*eta(1);            	% objective --> subjective map
s = a;                              % by definition

% prior (subjective)
m0 = mean(m(:,1));
if size(m,1) == 1
    s0 = 100000;
else
    s0 = (max(m)-min(m))/2.5;       % visual approx using Shi et al. (2013)    
end

% posteriors (subjective)
mh = (s^2*m0+s0^2*m)/(s^2+s0^2);
sh = sqrt((1/s^2+1/s0^2)^(-1));

% subjective --> objective map
muD = exp(m/eta(2))-1;             	% likelihood means
sigD = s*(exp(mu/eta(2))-1)/eta(2);	% likelihood SDs
muh = exp(mh/eta(2))-1;            	% posterior means
sigh = sh*(exp(mh/eta(2))-1)/eta(2);% posterior SDs
mu0 = exp(m0/eta(2))-1;            	% prior mean
sig0 = s0*(exp(m0/eta(2))-1)/eta(2);% prior SD

% Note: Because of logarithmic mapping, not necessary that lh > l or l0 
% in objective space (only in subjective space).

% distributions
prior = normpdf(t,mu0,sig0)';
[like,post] = deal(nan(length(t),length(mu)));
for e = 1:length(mu)
    like(:,e) = normpdf(t,muD(e),sigD(e));
    post(:,e) = normpdf(t,muh(e),sigh(e));
end

end