% matlab script to convert batch-analyzed qblade multi-reynolds cl cd 
% data into turbinesFoam usable profile files
% run the script in the same directory as the output files
% qblade output format should be aerodyn (nrel)

clear all; close all; clc;
list=dir('*.dat'); % get file names
list=struct2cell(list);
list=list(1,:)'; % get file names only as strings
nFiles = length(list); % get number of files
allData = cell(nFiles,4); % initialize cell that will store everything
AoA = transpose(-25:1:45); % desired AoA range
mkdir matlabOutput

for i=1:nFiles
    loc_ = strfind(list{i},'_');
    profileNameTemp=extractBetween(list{i},1,loc_(1)-1);
    allData{i,1}=string(profileNameTemp);
    profileReTemp=extractBetween(list{i},loc_(2)+3,loc_(3)-1);
    allData{i,2}=str2double (string(profileReTemp) )*1e6;
    clcdDataTemp = importdata(list{i},' ',14)           ;
    clcdDataTemp2 = clcdDataTemp.data    ;
%     clcdDataTemp2 = unique(clcdDataTemp.data,'rows')    ;
    allData{i,3} = interp1(clcdDataTemp2(:,1), clcdDataTemp2(:,2),AoA,'linear');
    allData{i,4} = interp1(clcdDataTemp2(:,1),clcdDataTemp2(:,3),AoA,'linear');
    
    allData{i,3} = fillmissing(allData{i,3},'linear');
    allData{i,4} = fillmissing(allData{i,4},'linear');
    
%     % plots for checking interpolation
%     % plot cl
%     figure; hold on;
%     plot(clcdDataTemp2(:,1),clcdDataTemp2(:,2))
%     plot(AoA,allData{i,3},'rx')
%     xlabel('AoA');ylabel('c_l');legend('qblade','interpolated')
%     % plot cd    
%     figure; hold on;
%     plot(clcdDataTemp2(:,1),clcdDataTemp2(:,3))
%     plot(AoA,allData{i,4},'rx')
%     xlabel('AoA');ylabel('c_d');legend('qblade','interpolated')
%     pause
%     close all
end

allData;
clear profileReTemp profileNameTemp loc_ i clcdDataTemp clcdDataTemp2
profileNames=unique( string(allData(:,1)) );
nProfiles=length(profileNames);
for i=1:nProfiles
    cellIndex = find([allData{:,1}] == profileNames(i));
    nRe = length(cellIndex);
    ReTemp = horzcat(allData{cellIndex,2});
    clTemp = horzcat(allData{cellIndex,3});
    cdTemp = horzcat(allData{cellIndex,4});
    
    [~,ReSort]=sort(ReTemp); % get order of ReTemp vector
    ReTemp=ReTemp(:,ReSort);
    clTemp=clTemp(:,ReSort);
    cdTemp=cdTemp(:,ReSort);
    
    fmtRe = [repmat('%11.4E   ',[1,size(ReTemp,2)]) '\n'];
    fmtReC = ['//( AOA         ' repmat('%11.4E   ',[1,size(ReTemp,2)]) ' )\n'];
    fmt =    ['( ' repmat('%11.4E   ',[1,size(clTemp,2)+1]) ' )\n'];

    fileRE = fopen(['matlabOutput/' profileNames{i} '_multiRe_Re.dat'],'w');
    fileCL = fopen(['matlabOutput/' profileNames{i} '_multiRe_cl.dat'],'w');
    fileCD = fopen(['matlabOutput/' profileNames{i} '_multiRe_cd.dat'],'w');

    fprintf(fileRE,['// ' datestr(now,'mm.dd.yyyy HH:MM') '\n']);
    fprintf(fileRE,fmtRe,ReTemp);

    fprintf(fileCL,['// ' datestr(now,'mm.dd.yyyy HH:MM') '\n']);
    fprintf(fileCL,fmtReC,ReTemp);
    fprintf(fileCL,fmt,[AoA,clTemp]');
    
    fprintf(fileCD,['// ' datestr(now,'mm.dd.yyyy HH:MM') '\n']);
    fprintf(fileCD,fmtReC,ReTemp);
    fprintf(fileCD,fmt,[AoA,cdTemp]');
    fclose('all');
end
