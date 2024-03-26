
function [Ap] = discreteeffortssharedtarget(Simhbdata,k,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill, hsb,fa,fm,target,bycatch)

%% Author --------------------------------------------------------------
% name: - Appilineni Kushal
% date - 7-25-2022
% email - akushal@ucdavis.edu

%% Assumptions
% 1. Catchability of all harvested fishes are same
% 2. Catchability of all bycatch species are same
% 3. Advantage and disadvantage of information is implicit i.e. both grow
% linearly with information exchange
% 4. Fishery opens only when bycatch> min cutoff bycatch and target species > min set by
% management
% 5. no recovery/seasonality
%inE = initial effort of all fisher, highest biomass target and bycatch are
% Here costs of fishing are high, costs are decided uniquely for each
% network
%emu = effective sensitivity of efforts to price

%simhbdata is conserved networks data after fishing free period
%simhbdata{1} = web;
%simhbdata{2} = fish;
%simhbdata{3} = biomass after fishing free steady state;
%simhbdata{4} = trophic level;

%% Description ---------------------------------------------------------
% 1.  Run the second stage of simulation after fishing free stage on the conserved foodwebs..
% 2. At the begining of fishing, a target and a bycatch is selected, TAC is calculated based on the biomasses from
% fishing free times. Target is selected highest biomass, bycatch is
% selected randomly
% 3. fishing occurs til the TAC is attained, both for harvest and bycatch
% 4. Furthermore, biomass of other species are monitored as bycatch. If biomass of a
% bycatch target falls below the cutoff of their biomass at t=0, policy implementation for bycatch is added (for
% conservation)
% 5. Bycatch policy - fishery closure until next season
%hill and hsb are given as input arrays
%baseline hsb and hill for species are 0.2 and 1.2
%hsb and hill inputs tell the values for trophic levels >4
%% SIMULATIONS -------------------------------------------------------------------------

    %% SETUP
   % N = 5; %no of fishers
   %k=randi(length(SimConsLin));
    web = Simhbdata{k}{1};
    
    %web=SimConsLin(k).topo.web; %randomly choose a conserved web from the list used in this study
    spe=length(web); % # of species? 
    fish = logical(Simhbdata{k}{2});
    B0 = Simhbdata{k}{3};
    %fish=SimConsLin(k).topo.fish;  %delineate the fish
    %ext=SimConsLin(k).free.ext; % delinate the fish that went extinct 
   % B0=SimConsLin(k).free.B; %biomass after the free stage (first 4000 time steps) 
    %qd=0;
    T = Simhbdata{k}{4};
    ihill = 1.2*(ones(spe,1));
    ihsb = 0.2*(ones(spe,1)); %initial values
    ihill(T>4) = hill;   %upper consumers have higher hill and hsb given by the inputs
    ihsb(T>4) = hsb;
    %T=SimConsLin(k).initial.Troph; %array of Trophic levels 
    %[r,K,y,e,Bs,c,h,ax_ar,Z,po,~]=setup_default(web,fish,hsb);
    [r,K,y,e,c,ax_ar,Z,po,~]=newsetup_default(web,fish);
    x=ax_ar.*(Z.^(-po.*(T-1))); %use the T to set up the metabolic rates
   % mu=0.3; %mu: market reactivity
    %co=0.01; %co: cost of fishing per unit of effort
    %cat=0.4; %catchability per unit of effort per unit of target biomass
    %cab=0.4; % catchability per unit of effort per unit of bycatch biomass
    %at=30;
    %bt=0.1;%0.01;
    %ab=-30;
    %bb= - 0.1;
   % cutofft = 0;  %cutoff for the permissible levels of target biomass depletion
    %cutoffb = 0;  %cutoff for the permissible levels of bycatch biomass depletion
    %price='linear';
   % pt=20;
   
    ui = false(spe,1);
    ui(target) = 1;
    target = ui;
    ui = false(spe,1);
    ui(bycatch) = 1;
    bycatch = ui;
    %%detectibility of n fishers.
    %de = sort(rand(1,N)); %detectibility is chosen as a random factor. when de = 1, then all effort goes towards target and if de = 0, all effort goes towards bycatch
    de = 1:N;
    de = 0.65 + sort(de/(4*N)); %increasing sequence of detectibility
    %% RUN

    tspan = [0,0.5,1]; %day to day integration
    %tspan=tspan/100;
    E0 = inE*ones(1,N); %set E0 for the simulation starting point as all fishers start with same effort = 1
    Tar0 = zeros(1,N*nharvest);%set initial harvest of individual target species by each fisher = 0
    By0 = zeros(1,N*nbycatch);%set initial bycatch of individual target species by each fisher = 0
    TS=[B0,E0,Tar0,By0];   %Biomass, effort, fishing season
    Btmin = B0(target)*cutofft; %minimum levels of biomass above which fishing is permitted
    options=odeset('RelTol',10^-8,'AbsTol',10^-8);

    Bcutoffeff = (costratio./de)*(sum(B0(target)));   %respective target biomass after which effort would go to 0 for fisher
   Bext = 10^(-6);
 
   
    for w=1:nseasons
        open = 1;       %fishing is open
        for days=1:365
            TS1 = TS;
            X0 = TS(size(TS,1),:);
            B0 = X0(1,1:spe);
           
            TAC = B0(target)-Btmin; %difference in current biomass and minimum levels of target
            
            if open == 1 && isempty(TAC(TAC<=0)) == 1  %none of the target species minimum is reached i.e. fishery open and none of bycatch is not reduced to 0
                [~,X] = ode45(@(t,X) d_seasonalcatchV2(t,X,N,de,x,y,r,K,e,ihill,c,ihsb,web,target,bycatch,cat,cab,Bext,fa,fm),tspan,X0,options);    %consecutive seasons
               
            else
                open = 0;
                E0 =zeros(1,N);     %effort goes to 0
                X0=[B0,E0,Tar0,By0];   %Biomass, effort, fishing season
                [~,X] = ode45(@(t,X) newseasonclosure(t,X,N,nharvest,nbycatch,x,y,r,K,e,ihill,c,ihsb,web,Bext,fa,fm),tspan,X0,options);    %consecutive seasons
            end
            B0=X(3,1:spe);
            E0=X(3,spe+1:spe+N);
            Tar0 = X(3,spe+1+N:spe+(nharvest+ 1)*N);
            By0 = X(3,spe+1+(nharvest+ 1)*N:spe+(nharvest + nbycatch+ 1)*N);
            Tar0(Tar0<Bext)=0;
            By0(By0<Bext)=0;
            B0(B0<Bext)=0;
           
            Echeck = find((B0(target) - Bcutoffeff)<=0);
            E0(Echeck) = 0;
            TAC = B0(target)-Btmin; %difference in current biomass and minimum levels of target
            if days==365 && isempty(TAC(TAC<=0)) == 1%last day of fishing season after closure
                Echeck = find((B0(target) - Bcutoffeff)>0); % biomass such that it's profitable
                E0(Echeck) = inE;
               
            end
            X0=[B0,E0,Tar0,By0];  %time series data of that 1 day
            TS = [TS1;X0];   
        end
    end

    Ap = cell(1,2);
    Ap{1} =TS;
    Stats=zeros(1,2);
    Stats(1,1)= find(target); %target index
    Stats(1,2)=find(bycatch); %bycatch index
    
    Ap{2} = Stats;
        


end