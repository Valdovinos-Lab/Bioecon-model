%% webgeneration.m

%% Author -------------------------------------------------------------
% name: Valentin Cocco
% mail: valentin.cocco@ens.fr
% creation: 4-2-2018

%% Modified by -------------------------------------------------------------
% name: Appilineni Kushal
% mail: akushal@ucdavis.edu 
% creation: 7-7-2023


%% Description --------------------------------------------------------
% Using the niche model, create a large number of foodwebs represented by their adjacency matrix. The foodwebs are saved
% in a MAT-file. The description of the foodweb also include initial biomass, fish identities.
% Calls:
%   - NicheModel
% Inputs:
%   - n: number of foodwebs produced.
%   - spe: number of species in the foodweb
%   - con: connectance of the foodweb
%   - err: authorized relative error of connectance
% Outputs:
%   - web: adjacency matrix of a foodweb
%   - fish: logical vector indicating fish identity or not
%   - B0: initial biomass


%%
n=3000; %number of networks to create
spe=30;
con=0.15;
err=0.025;

webs=cell(n,3);

for i=1:n
    sprintf('Web %d/%d', i, n)
    
    web=NicheModel(spe,con,err);
    
    % Trophic level T calculated with the algebric method of Levine, 1980.
    W=zeros(spe);
    nocanweb=web;
    nocanweb(logical(eye(spe)))=0;
    prey=sum(nocanweb,2)*ones(1,spe);
    W(prey~=0)=nocanweb(prey~=0)./prey(prey~=0); %W_i_j: 1/number of preys of i if i eats j; 0 otherwise
    T=(inv(eye(spe)-W))*ones(spe,1);
        
    fish=false(spe,1);
    fishp=0.6;
    bernoulli=rand(spe,1);
    fish(T>=3 & bernoulli<fishp)=true;

    B0=0.05+0.95*rand(spe,1);
    
    webs{i,1}=web;   %stores the food web
    webs{i,2}=fish;  %stores the fish species
    webs{i,3}=B0;    %stores the initial biomasses
    
%     % Save in a .csv file
%     tableweb=array2table([web,fish,B0]);
%     cd('Data/WebStock')
%     filename=sprintf('Web%03d.csv',i);
%     writetable(tableweb,filename,'WriteVariableNames',false,'Delimiter',',')
%     cd('../..')
end

% Save in a MAT-file
%cd('Data')
save('Webs3000.mat','webs') %Change name of file to create new webs
