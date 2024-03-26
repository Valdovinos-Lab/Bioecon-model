 %% newwebfishingfree.m
clear
%% Previous Author --------------------------------------------------------------
% when: 8-14-2019
% who: Paul Glaum (prglaum@umich.edu)
% what: used to run simulations without fishing effort for 500 different foodweb structures.

%% Description ---------------------------------------------------------
% Using the 3000-foodweb-structure stock, run simulations for 4000 timesteps without fishing effort locally.
% This serves as an example of the ecological dynamics of the aquatic
% networks without fishing. 

% Calls:
%   - webproperties.m
%   - newsetup_default.m
%   - newdifferentialresp.m

%% Updates -----------------------------------------------------------------
% when: 7-7-2023
% who: Appilineni Kushal (akushal@ucdavis.edu)
% what: include maintenence respiration losses and 
% assimilated carbon for production from Nadja Kath paper
% also, set half saturation biomass and hill coefficient in this function directly
%% 1. UPLOADING OF THE NETWORK STRUCTURE
load('Webs3000.mat')
%con=0.15; %(Berlow2009)
%err=0.025; %(Williams2000: 0.03)
nSim=1; %number of simulations
spe=30; %number of species
t1=4000; %last timestep
tspan=0:t1;
Fail = cell(1,3000);   %storing the failed cases along with their index

newUT4Fa4Fm4H12Hb2data = cell(1,3000);
parfor k=1:3000
    tot=zeros(spe,1); 
   % sprintf('Simulation %d/%d',k,nSim)
    hill = 1.2*ones(spe,1);
    hsb = 0.2*ones(spe,1);
    web=webs{k,1};
    fish=webs{k,2};
    B0=webs{k,3};
    basal=find(sum(web,2)==0);    %%%%%shows species with no prey
    %% 2. CALCULATION OF THE INITIAL STRUCTURAL PROPERTIES
    [tmp, T]=webproperties(web);
    Fish=nnz(fish)/spe;
    properties=[tmp,Fish];
   hill(T>4) = 1.2;   %upper consumers can have higher hill and hsb
    hsb(T>4) = 0.2;
    %% 3. SET THE BIOLOGICAL PARAMETERS
    [r,K,y,e,c,ax_ar,Z,po,Bext]=newsetup_default(web,fish);
    x=ax_ar.*(Z.^(-po.*(T-1)));
    fa = 0.4;
    fm = 0.4;
    %% 4. ECONOMIC PARAMETERS (These are irrelevant for fishing free simulations)
    mu=0;
    co=1;
    ca=0;
    a=1;
    b=0;
    price='linear'; 

    %% 5. ATN MODEL
    E0=[];  %efforts are 0, thus economic dynamics are irrelevant
    harv=false(spe,1);
    X0=[B0,E0];
    options=odeset('RelTol',10^-8,'AbsTol',10^-8);
   % [t,X] = ode45(@(t,X) newdifferential(t,X,x,y,r,K,e,hill,c,hsb,web,harv,mu,co,ca,a,b,price,Bext),tspan,X0,options);
    [t,X] = ode45(@(t,X) newdifferentialresp(t,X,x,y,r,K,e,hill,c,hsb,web,harv,mu,co,ca,a,b,price,Bext,fa,fm),tspan,X0,options);
    B=X(:,1:spe);
    E=X(:,spe+1:end);
    B(B<Bext)=0;
    E(E<0)=0;
    X=[B,E];
    
    %% no of species > 20
    finalbiomass = X(4001,:);
    ok = 1;
    if nnz(finalbiomass)>=20 && nnz(finalbiomass(fish))>0  %contains 20 species atleast and atleast one fish present
        extinct = find(finalbiomass == 0);
        df = true(1,spe);
        df(extinct) = 0; %array containing the index of extinct species
        for o = 1:length(extinct) %removing extinct species from adjacency matrix by setting their connections to 0
            web(extinct(1,o),:) = 0;
            web(:,extinct(1,o)) = 0;
        end
        predator=sum(web,1)'; %number of predator per species
        prey=sum(web,2); %number of prey per species
        link=predator+prey; %number of links per species
        %% isolated species
        if all(link(df))==0  
            ok=0;
            Fail{k} = [Fail{k} 1];
        else
            %% 2. Is there any species not connected to a basal species ? 
            %tot=zeros(spe,1); %tot(i)=1 : species i connected to a basal species
            newbasal = basal.*(X(4001,basal)>0)';
            basal  = newbasal(newbasal>0); %removing any basal species that might have gone extinct
            tmp=basal; %tmp: last species connected to basal species found
            tmpweb=web;
            tot(tmp)=1;
            while isempty(tmp)==0 %when tmp=[], there is no more species connected to basal ones to find
                [i,~]=find(tmpweb(:,tmp)~=0);
                tmpweb(:,tmp)=0;
                tmpweb(tmp,:)=0; %the last species "disappear" from the web to avoid infinite loops (i eats j, that eats k, that eats i...)
                tmp=unique(i); %i predators of tmp --> i are connected to basal species
                tot(tmp)=1;       
            end
            %% 

            if all(tot(df))==0 %at least one species is not connected to a basal
                ok=0;
                 Fail{k} = [Fail{k} 2]
            else
                %% 3. Is there more than one independant network ?
                allZeroIndices = all(web == 0, 1) & all(web == 0, 2)';
                connectedweb = web(~allZeroIndices, ~allZeroIndices);
                z=isConnected(connectedweb);
                %% 
                if z==0
                    ok=0;
                     Fail{k} = [Fail{k} 3];
%                 else
%                 %% 4. Is the actual connectance close enough from con ?
%                     con_web=nnz(connectedweb)/(nnz(finalbiomass))^2;
%                     con_min=(1-err)*con;
%                     con_max=(1+err)*con;
%                 %% 
%                     if con_web<con_min || con_web>con_max
%                         ok=0;
%                         Fail{k} = [Fail{k} 4];
%                     end 
               end
            end
        end
    else
        ok =0;  
    end
    %create a list of things that need to be stored
    %web, fish species, initial biomass, trohpic level array
    re = false(spe,1);
    re(fish) = 1;
    Ap = cell(1,5);
   
    if ok ~= 0
       Ap{1} = web; %store the web
       Ap{2} = re.*(X(4001,:)>0)'; %store the fish
       Ap{3} = X(4001,:); %biomass
       Ap{4} = T;  
       Ap{5} = X
    else
        Ap{2} = fish;
        Ap{3} =  X(4001,:); %biomass
        Ap{4} = T;
        %Ap{1} = tot;
    end
    newUT4Fa4Fm4H12Hb2data{k} = Ap;
end
save newUT4Fa4Fm4H12Hb2data.mat newUT4Fa4Fm4H17Hb5data -v7.3
