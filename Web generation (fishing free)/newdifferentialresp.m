%% ATNmodel.m

%% Previous authors -------------------------------------------------------
% Valentin Cocco

%% Last update ------------------------------------------------------------
% when: 7-7-2023
% who: Appilineni Kushal
% mail: akushal@ucdavis.edu
% added respiration loses and assimilated carbon used for production based
% on Nadja Kath 2018 work

%% Description ------------------------------------------------------------
% Determine the temporal derivative of B and E
% Called by:
%   - newwebfishingfree.m
% Inputs:
%   - X: [B;E]
%   - x: array of metabolic rates
%   - y: array of maximum assimilation rates
%   - r: array of growth rates
%   - c: array of intraspecific predator interference
%   - K: global carrying capacity
%   - e: matrix of metabolization efficiency
%   - Bs: array of half-saturation density
%   - nicheweb: 
%   - harvest: logical array. indicates which trophic species are harvested
%   - mu: market openness
%   - co: cost of fishing per unit of effort
%   - ca: catchability per unit of effort per unit of biomass
%   - a, b: price parameters
%   - price: indicate the price model chosen : 'linear', 'isoelastic', 'nl-ni'
%   - Bext: extinctionthreshold
%   - fa: Fraction of assimilated carbonused for production
%   - fm: Factor for maintenence respiration

%%
function [dXdt]=newdifferentialresp(t,X,x,y,r,K,e,h,c,Bs,nicheweb,harvest,mu,co,ca,a,b,price,Bext,fa,fm)
spe=length(nicheweb);
B=X(1:spe);
E=X(spe+1:end);
B(B<Bext)=0;
E(E<0)=0;

%% ----------------------------------------------------------------------------
% BIOENERGETIC MODEL ----------------------------------------------------------
% -----------------------------------------------------------------------------

%% 1. NPP: Net Primary Production
basal=(sum(nicheweb,2)==0);
NPP=r.*(1-sum(B(basal))./K).*B;

%% 2. MetaLoss: Metabolic Loss (mortality, energy exependiture for basal metabolism, activity, thermoregulation...)
MetaLoss=(fm*(~r)).*x.*B; %metabolic losses for everyone but basal

%% 3. F: Functional Response
w=zeros(spe,1);
nonbasal=~basal;
w(nonbasal)=1./sum(nicheweb(nonbasal,B~=0),2); %only the remaining species are taken into account
w(w==Inf)=0; %predators whose all preys have a null biomass
w=w*ones(1,spe);
w=w.*nicheweb; %w_i_j: preference of i for j amongst its preys

wBh=w.*(ones(spe,1)*B').^(ones(spe,1)*h');
Bsh=(ones(spe,1)*Bs').^(ones(spe,1)*h');
cBBsh=(c.*B*ones(1,spe)).*Bsh;
sumwBh=sum(wBh,2)*ones(1,spe);
F=wBh./(Bsh+cBBsh+sumwBh); %Bsh is always different from zero --> no problem of division

%% 4. TrophLoss and TrophGain: Loss and Gain of organic matter by consumer-resource interaction
gain=(x.*y.*B)*ones(1,spe).*F;
gain1=((fa*(~r)).*x.*y.*B)*ones(1,spe).*F; %including respiration in gains
loss=gain./e;
TrophGain=sum(gain1,2);
TrophLoss=sum(loss,1)';

%% 5. HarvLoss: Loss by fishing
HarvLoss=zeros(spe,1);
HarvLoss(harvest)=ca*E.*B(harvest);

%% 6. dBdt: array of temporel derivatives of B_i
dBdt=NPP - MetaLoss - TrophLoss + TrophGain - HarvLoss;

%% -----------------------------------------------------------------------------
% ECONOMIC MODEL ---------------------------------------------------------------
% ------------------------------------------------------------------------------
%% 7. p: price
Y=HarvLoss(harvest);
if price=='linear'
    p=a*(1-b*Y);
elseif price=='isoelastic'
    p=a/Y^b;
    p(p==Inf)=0;
else %non linear non isoelastic: 'nl-ni'
    p=a/(1+b*Y);
    p(p==Inf)=0;
end

%% 8. dEdt
dEdt=mu*(ca*p.*B(harvest).*E-co*E);


%% --------------------------------------------------------------------------------
dXdt=[dBdt;dEdt];
