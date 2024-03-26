%% newsetup_default.m

%% Previous authors -------------------------------------------------------
% when: 4-2-2018
% who: Valentin Cocco
% mail: valentin.cocco@ens.fr
% what: personal writing for constant effort framework

%% Last update ------------------------------------------------------------
% when: 7-7-2023
% who: Appilineni Kushal (akushal@ucdavis.edu)
% what: doesn't return Hill's coefficient or Half saturation biomass, so
% that we can set it according to our convenience in the
% newwebfishingfree.m 

%% Description ------------------------------------------------------------
% Set the biodynamical DEFAULT parameters of the foodweb.
% Used for  the constant effort framework
% Inputs:
%   - nicheweb: adjacency matrix (row i eats column j)
%   - fish: logical array; indicate which species are fishes
% Called by:
%   - newwebfishingfree.m

%%
function [r,K,y,e,c,ax_ar,Z,po,Bext]=newsetup_default(nicheweb,fish)

%% ------------------------------------------------------------------------
% 1. DYNAMICAL PARAMETERS -------------------------------------------------
% -------------------------------------------------------------------------
% r: growth rate of i (non zero for basal species only)
% K: global carrying capacity of i (for basal species only)
% y: maximum consumption rate of i eating j
% e: assimilation efficiency of i eating j
% Bs: half saturation biomass for i eating j
% c: interference of predator i eating j with other predators of j
% h: Hill coefficient

%% r: growth rate ---------------------------------------------------------------------------
% constant r=1 for every basal species
spe=length(nicheweb);
r=zeros(spe,1);
basal=(sum(nicheweb,2)==0);
r(basal)=1;

%% K: global carrying capacity ----------------------------------------------------------------
K=ones(spe,1);

%% x: metabolic rate --------------------------------------------------------------------------
nonbasal=(sum(nicheweb,2)~=0);

ax_ar=zeros(spe,1);
ax_ar(nonbasal)=0.314;
ax_ar(fish)=0.88;

Z=100;

po=0.25;

%% y: maximal consumption rate per unit or metabolic rate -------------------------------------
y=zeros(spe,1);

%y(nonbasal)=10;
y(nonbasal) = 8;
y(fish) = 4;

%% e: assimilation efficiency -----------------------------------------------------------------
% carnivorous/herbivorous
%i eats non basal (carnivorous)
e=0.85*ones(spe);
%i eats basal (herbivorous)
e(:,basal)=0.45;

%% Bs: half saturation biomass ---------------------------------------------------------------
%Bs=hsb*ones(spe,1);

%% c: feeding interference coefficient -------------------------------------------------------
% Holling-type II or III
c=zeros(spe,1);

%% h: Hill coefficient ------------------------------------------------------------------------
%h=hill;

%% Bext: extinction threshold -----------------------------------------------------------------
Bext=10^-6;
