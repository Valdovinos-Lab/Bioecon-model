%% ATNmodel.m

%% Previous authors -------------------------------------------------------
% Program by: Rosalyn Rael
% Modified by Barbara Bauer

%% Last update ------------------------------------------------------------
% when: 7-22-2022
% who: Appilineni Kushal
% mail: akushal@ucdavis.edu
% what: personal writing

%% Description ------------------------------------------------------------
% Determine the temporal derivative of B and E
% Called by:
%   - webdriver.m
% Inputs:
%   - X: [B;E]
%   - Btmin: Minimum biomasses (target) for fishing to be open
%   - de: detectibility/information aquired by the n fishers
%   - x: array of metabolic rates
%   - y: array of maximum assimilation rates
%   - r: array of growth rates
%   - c: array of intraspecific predator interference
%   - K: global carrying capacity
%   - e: matrix of metabolization efficiency
%   - Bs: array of half-saturation density
%   - nicheweb: 
%   - target: logical array. indicates which trophic species are targets
%   - bycatch: logical array. indicates which trophic species are bycatch
%   - mu: market openness
%   - co: cost of fishing per unit of effort
%   - cat: catchability per unit of effort per unit of target biomass
%   - cab: catchability per unit of effort per unit of bycatch biomass
%   - at, bt: price parameters for target
%   - ab, bb: price parameters for bycatch
%   - price: assumed constant
%   - Bext: extinctionthreshold

%added hill and hsb for higher trophic levels
%added respiration loses
%%
function [dXdt]=newseasonclosure(t,X,N,nharvest,nbycatch,x,y,r,K,e,h,c,Bs,nicheweb,Bext,fa,fm)
spe=length(nicheweb);
B=X(1:spe);
E=X(spe+1:spe+N);
Tar = X(spe+1+N:spe + (nharvest+ 1)*N);
By = X(spe+1+(nharvest+ 1)*N:spe+ (nharvest + nbycatch+ 1)*N);
B(B<Bext)=0;
Tar(Tar<Bext)=0; 
By(By<Bext)=0;
E(E<0)=0;
%E(E>10000)=10000;
%open = X(end);  %an indicator that fishery is open for the year or not
%% ----------------------------------------------------------------------------
% BIOENERGETIC MODEL ----------------------------------------------------------
% -----------------------------------------------------------------------------

%% 1. NPP: Net Primary Production
basal=(sum(nicheweb,2)==0);
NPP=r.*(1-sum(B(basal))./K).*B;

%% 2. MetaLoss: Metabolic Loss (mortality, energy exependiture for basal metabolism, activity, thermoregulation...)
MetaLoss=(fm*(~r)).*x.*B; 

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


%% 6. dBdt: array of temporel derivatives of B_i
dBdt=NPP - MetaLoss - TrophLoss + TrophGain;
dTardt = zeros(N*nharvest,1);
dBydt = zeros(N*nbycatch,1);

%% -----------------------------------------------------------------------------
% ECONOMIC MODEL ---------------------------------------------------------------
% ------------------------------------------------------------------------------

dEdt =zeros(N,1);
dXdt=[dBdt;dEdt;dTardt;dBydt];





