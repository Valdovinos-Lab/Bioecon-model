%% ATNmodel.m

%% Previous authors -------------------------------------------------------
% Program by: Rosalyn Rael
% Modified by Paul Glaum

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
%   - emu: effective market sensitivity
%   - co: cost of fishing per unit of effort
%   - cat: catchability per unit of effort per unit of target biomass
%   - cab: catchability per unit of effort per unit of bycatch biomass
%   - at, bt: price parameters for target
%   - ab, bb: price parameters for bycatch
%   - costratio: effective fishing costs/profits at t=0
%   - Bext: extinctionthreshold
% introduced hill and hsb for higher trophic levels
%added respiration loses and maintence 

%%
function [dXdt]=d_seasonalcatchV2(t,X,N,de,x,y,r,K,e,h,c,Bs,web,target,bycatch,cat, cab,Bext,fa,fm)
spe=length(web);
nharvest = sum(target);
nbycatch = sum(bycatch);
B=X(1:spe);
E=X(spe+1:spe+N);
%Tar = X(spe+1+N:spe+2*N);
%By = X(spe+1+2*N:spe+3*N);
Tar = X(spe+1+N:spe + (nharvest+ 1)*N);
By = X(spe+1+(nharvest+ 1)*N:spe+ (nharvest + nbycatch+ 1)*N);
B(B<Bext)=0;
Tar(Tar<Bext)=0; 
By(By<Bext)=0;
E(E<10^(-8))=0;
%E(E>inE) = inE;
%E(E>10000)=10000;
%open = X(end);  %an indicator that fishery is open for the year or not
%% ----------------------------------------------------------------------------
% BIOENERGETIC MODEL ----------------------------------------------------------
% -----------------------------------------------------------------------------

%% 1. NPP: Net Primary Production
basal=(sum(web,2)==0);
NPP=r.*(1-sum(B(basal))./K).*B;

%% 2. MetaLoss: Metabolic Loss (mortality, energy exependiture for basal metabolism, activity, thermoregulation...)
MetaLoss=(fm*(~r)).*x.*B;

%% 3. F: Functional Response
w=zeros(spe,1);
nonbasal=~basal;
w(nonbasal)=1./sum(web(nonbasal,B~=0),2); %only the remaining species are taken into account
w(w==Inf)=0; %predators whose all preys have a null biomass
w=w*ones(1,spe);
w=w.*web; %w_i_j: preference of i for j amongst its preys

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
%TAC = B(target) - Btmin'; %difference in current biomass and minimum levels of target
%if isempty(TAC(TAC<=0)) == 1 && isempty(bg(bg==0))==1 %none of the target species minimum is reached i.e. fishery open and none of bycatch is not reduced to 0
Totalefforttarget = sum(de.*E'); %sum of all the efforts by N fishers towards target
Totaleffortbycatch = sum((1-de).*E'); %sum of all the efforts by N fishers towards bycatch
HarvLoss(target)=cat*Totalefforttarget.*B(target);
HarvLoss(bycatch)=cab*Totaleffortbycatch.*B(bycatch);
%end

%% 6. dBdt: array of temporel derivatives of B_i
dBdt=NPP - MetaLoss - TrophLoss + TrophGain - HarvLoss;
%dTardt = (cat*(de.*E')*sum(B(target)))';
%dBydt = (cab*((1-de).*E')*sum(B(bycatch)))';
Tbiomass = B(target);
Bybiomass = B(bycatch);
Tch=0;
for i=1:nharvest         %change of individual target by individual fishers
    Ph = (cat*(de.*E')*Tbiomass(i))';
    Tch = [Tch;Ph];
end
Tch(1) = [];
dTardt = Tch;
Tch=0;
for i=1:nbycatch        %change of individual bycatch by individual fishers
    Ph = (cab*((1-de).*E')*Bybiomass(i))';
    Tch = [Tch;Ph];
end
Tch(1) = [];
dBydt = Tch;
%% -----------------------------------------------------------------------------
% ECONOMIC MODEL ---------------------------------------------------------------
% ------------------------------------------------------------------------------
%if isempty(TAC(TAC<=0)) == 1 && isempty(bg(bg==0))==1 %if fishery is open
%% 7. pt: price of target and pb = price of bycatch
    %Ytarget=sum(HarvLoss(target));  %total yield of target species
   % Ybycatch=sum(HarvLoss(bycatch));%%total yield of bycatch species 
   % if price=='linear'
     %   pt=at*(1-bt*Ytarget);  %price equation for target
     %   pb=ab*(1-bb*Ybycatch);  %price equation for bycatch
   % elseif price=='isoelastic'
    %    p=a/Y^b;
    %    p(p==Inf)=0;
   % else %non linear non isoelastic: 'nl-ni'
    %    p=a/(1+b*Y);
    %    p(p==Inf)=0;
  %  end

%% 8. dEdt
%result = E>=inE;
%dEdt=emu*(E).*(cat*(de*sum(B(target))) - effectivecosts - qd*E')';  %cost free information
%dEdt(result) = min(0,dEdt(result));
dEdt = zeros(length(E),1);

%% --------------------------------------------------------------------------------
dXdt=[dBdt;dEdt;dTardt;dBydt];





