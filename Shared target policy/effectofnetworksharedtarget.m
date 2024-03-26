function [] = effectofnetworksharedtarget(start)
    load('newwebsUT4Fa4Fm4H12Hb2data.mat', 'A')
    load('dataforrandomtarget.mat', 'dataforrandomtarget')
    N=5;
    %qd=0;
    %emu=0.3*20; %mu: market reactivity
  % emu = 1;
 % emu = 50;
   % co=0.01; %co: cost of fishing per unit of effort
    costratio = 0.25;
    cat=0.002; %catchability per unit of effort per unit of target biomass
    cab=0.002;
    cutofft=[0.15 0.3 0.45 0.6];
    cutoffb = 0;
   % pt=20;
    nharvest=1;
    nbycatch=1;
    inE=1;
   hsb = 0.2;  %half saturation biomass
   fa = 0.4;
   fm = 0.4;
   hill = 1.2;
    nseasons = 5;
    ty= start(1);
    ty2=ty + 39;
    justtheparforpartsharedtarget(ty,ty2,A,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,dataforrandomtarget);
    
end
