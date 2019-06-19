%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code developed by Mahshid Ghanbari (mahshid.ghanbari@colostate.edu)

% Please reference to: Ghanbari, M., Arabi, M., Obeysekera,J., & Sweet, W. (2019.
% A coherentstatistical model for coastalfloodfrequency analysis under nonstationary
% sea level conditions.Earth's Future,7,162–177.
% https://doi.org/10.1029/2018EF001089

% Please contact Mahshid Ghanbari (mahshid.ghanbari@colostate.edu) with any issue.
%  Inputs:
%  Maximum and mean Daily Water level data ( The output from Daily_WaterLevel_Data.m)
%  Major and Minor flood threshold (https://tidesandcurrents.noaa.gov/publications/techrpt86_PaP_of_HTFlooding.pdf)
%  Station's Name

%  Outputs:
%  (1)Historical and future flood return level interval curves
%     Table: Return_Period_Level. The first column is annual return period and the rest of columns are 
%     retrun levels under 0 to 2 ft sea level rise (e.g., 0 ft, 0.1 ft, 0.2 ft, ...)
%     The figure will be seved in Outputs folder and insert the station name
%  (2)Future versus current flood return period
%     Table: Future_Return_Period. The first column is selected return
%     perido and the rest of columns are future return period under 
%     0 to 2 ft sea level rise (e.g., 0 ft, 0.1 ft, 0.2 ft, ...)
%  (3)Historical and future frequncy of minor and major flooding
%     Table: Minor_Major_Flood. The retrun period of major flooding and
%     annual frequncy of minor flooding under 0 to 2 ft sea level rise
%     (e.g., 0 ft, 0.1 ft, 0.2 ft, ...). 


% HOW TO RUN THE Code :
% (1) Place hourly observed sea level data(station name.csv) in the main
%     floder (e.g., Battery.csv)
% (2) Run Daily_WaterLevel_Data.m (The mean and maximum daily values
%      will be saved in Inputs folder. e.g.,Battery_DailyWaterLevel.xlsx )
% (3) run Run_Mixture_Model.m after inserting Minor and major flood threshold
%     and the desired return period

% Station_name={'Boston', 'Battery','GrandIsle', 'SanFrancisco'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;
clearvars -except Station_name

if exist('Station_name','var')
inputs.Station_name= Station_name;
else
  inputs.Station_name=input('Insert the station name :','s');  
end

%% User Input
inputs.Minor_Th= 1.84; % Insert minor flood threshold [ft]
inputs.Major_Th= 4.03; % Insert minor flood threshold [ft]

inputs.Return_Period=[1,50,100,200,500,1000]; % Insert the desired return period

%% Reading the mean and maximum daily values from Inputs folder
if ispc
    cd([pwd,'\Inputs'])
else
    cd([pwd,'/Inputs'])
end

inputs.DayMax_byYear=xlsread(strcat(inputs.Station_name,'_DailyWaterLevel.xlsx'),'Max');
inputs.DayMean_byYear=xlsread(strcat(inputs.Station_name,'_DailyWaterLevel.xlsx'),'Mean');
cd ..

[Minor_Major_Flood,Return_Period_Level,Future_Return_Period]=Mixture_Model(inputs);
clear Station_name inputs