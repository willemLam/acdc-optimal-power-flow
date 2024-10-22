
function [PV_perun, PV_solarmax, Pac_facade] = readFromMeteoBox()

%%
    reading_port =39019;
    [data] = readFromResources(reading_port, 2^12 );
    timestamp = data.Info.Timestamp;
    GHI = data.Data.GHI;
    GI = data.Data.GI;
    GVI = data.Data.GVI;
    Tair = 10;%data.Data.Tair;
    MatlabTime = epochToMatlabTime(timestamp);
    UTC = 2;

    [PV_perun, PV_solarmax,PV_facade] = PV_production(GHI, GI, GVI, Tair, MatlabTime, UTC);
   
    PV_perun
    PV_facade+PV_solarmax
  %%  
    PV_perun(isnan(PV_perun)) = 0 ;
    PV_solarmax(isnan(PV_solarmax)) = 0 ;
    Pac_facade(isnan(Pac_facade)) = 0 ;


end