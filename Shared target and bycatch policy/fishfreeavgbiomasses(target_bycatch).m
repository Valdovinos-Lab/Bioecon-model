load("")
load('dataforrandomtarget.mat')
t = 4001;
Fishfreeavgbiomass = cell(1,400); %calculate fish free target and bycatch biomass avg over last year
parfor i=1:400
    ui = dataforrandomtarget{i}; %target bycatch combinations for that network
    op = zeros(length(ui),2); %containing average target bycatch values
    Ts = A{i}{5}; %time series
    for i1=1:length(ui)
        target = ui(i1,1);
        bycatch = ui(i1,2);
        op(ui,:) = [mean(Ts(t-364:t,target))  mean(t-364:t,bycatch)]; 
    end
    Fishfreeavgbiomass{i} = op;
end