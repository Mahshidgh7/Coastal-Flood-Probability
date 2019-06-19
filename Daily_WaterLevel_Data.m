clear all; close all; clc;
%% User Input
Station_name= 'Battery'; % Insert the station name
%%
if ispc
    cd([pwd,'\Inputs'])
else
    cd([pwd,'/Inputs'])
end
filename = strcat(Station_name,'.csv');
fid=fopen(filename);
data = textscan(fid, '%s %s', 'Delimiter',',');
Time=data{1,1};
formatIn = 'dd-mm-yyyy ';
final_Data(:,1)= str2double(string(data{1,2} (:,1)));
final_Data(:,2)=datenum(Time(:,1),formatIn);
ctr=0;
for i = min(final_Data(:,2)): max(final_Data(:,2))
    ctr=ctr+1;
    index=find(final_Data(:,2)==i);
    if length(index)>10
        Daily_height(ctr,3)=nanmax(final_Data(index,1));
        Daily_height(ctr,4)=nanmean(final_Data(index,1));
    else
        Daily_height(ctr,3)=NaN;
        Daily_height(ctr,4)=NaN;
    end
    Daily_height(ctr,1)=i;
    Daily_height(ctr,2)=str2double(datestr(i,'yyyy'));
end
ctr=0;
for i = min(Daily_height(:,2)):max(Daily_height(:,2))
    clear Max_Daily_Year1 Mean_Daily_Year1
    ctr=ctr+1;
index1=find(Daily_height(:,2)==i);
Max_Daily_Year1=Daily_height(index1,3);
Mean_Daily_Year1=Daily_height(index1,4);
Max_Daily_Year1(isnan(Max_Daily_Year1))=[];
Mean_Daily_Year1(isnan(Mean_Daily_Year1))=[];
if length(Max_Daily_Year1)~=366 && length(Max_Daily_Year1)>1
    index2=(max(find(Max_Daily_Year1==Max_Daily_Year1(end))))+1;
    Max_Daily_Year1(index2:366)=NaN;
elseif length(Max_Daily_Year1)< 1
   Max_Daily_Year1(1:366)=NaN;
end
if length(Mean_Daily_Year1)~=366 && length(Mean_Daily_Year1)>1
    index3=(max(find(Mean_Daily_Year1==Mean_Daily_Year1(end))))+1;
    Mean_Daily_Year1(index3:366)=NaN;
elseif length(Mean_Daily_Year1)< 1
    Mean_Daily_Year1(1:366)=NaN;
end
Max_Daily_Year(:,ctr)=Max_Daily_Year1;
Mean_Daily_Year(:,ctr)=Mean_Daily_Year1;
end

if ~exist('Inputs','dir')
    mkdir('Inputs')
end

xlswrite(strcat(Station_name,'_DailyWaterLevel.xlsx'),Max_Daily_Year,'Max')
xlswrite(strcat(Station_name,'_DailyWaterLevel.xlsx'),Mean_Daily_Year,'Mean')

cd ..