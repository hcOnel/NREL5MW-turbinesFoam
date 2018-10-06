% matlab script to convert batch-analyzed qblade multi-reynolds cl cd 
% data into turbinesFoam usable profile files
% run the script in the same directory as the output files
% qblade output format should be aerodyn (nrel)
% v4: octave compatible

clear all; close all; clc;

willPlot=input("Want to see interpolation on plots? 1 or 0 \n");
if ~(willPlot==1|willPlot==0)
	error("Invalid input. Terminating.");
end

qBladePolarSubdir = "qBladePolarOutput/";
list=dir([qBladePolarSubdir '*.dat']); % get file names
list=struct2cell(list);
list=list(1,:)'; % get file names only as strings
nFiles = length(list); % get number of files
allData = cell(nFiles,4); % initialize cell that will store everything
AoA = transpose(-25:1:45); % desired AoA range
mkdir convertedClCdData

for i=1:nFiles
    loc_ = strfind(list{i},'_');
    profileNameTemp=list{i}(1:loc_(1)-1);
    allData{i,1}=profileNameTemp;
    profileReTemp=list{i}(loc_(2)+3:loc_(3)-1);
    allData{i,2}=str2num(profileReTemp)*1e6;
    clcdDataTemp = importdata([qBladePolarSubdir list{i}],' ',14)           ;
    clcdDataTemp2 = clcdDataTemp.data    ;
%     clcdDataTemp2 = unique(clcdDataTemp.data,'rows')    ;
    allData{i,3} = interp1(clcdDataTemp2(:,1), clcdDataTemp2(:,2),AoA,'linear','extrap');
    allData{i,4} = interp1(clcdDataTemp2(:,1),clcdDataTemp2(:,3),AoA,'linear','extrap');
    
	if willPlot==1
		% plots for checking interpolation
		% plot cl
		figure; hold on;
		plot(clcdDataTemp2(:,1),clcdDataTemp2(:,2))
		plot(AoA,allData{i,3},'rx')
		xlabel('AoA');ylabel('c_l');%legend('qblade','interpolated')
		title(profileNameTemp)
		pause(0.5)
		% plot cd    
		figure; hold on;
		plot(clcdDataTemp2(:,1),clcdDataTemp2(:,3))
		plot(AoA,allData{i,4},'rx')
		xlabel('AoA');ylabel('c_d');%legend('qblade','interpolated')
		title(profileNameTemp)
		pause(0.5)
		close all
	end
end

allData;
clear profileReTemp profileNameTemp loc_ i clcdDataTemp clcdDataTemp2
profileNames=unique( (allData(:,1)) )
nProfiles=length(profileNames);
for i=1:nProfiles

	% find indices of each profile
    cellIndex = strfind(allData(:,1),profileNames{i});
    for j=1:length(cellIndex)
		if isempty(cellIndex{j})==1
			cellIndex{j}=0;
		end
    end
    cellIndex = cell2mat(cellIndex);
    cellIndex=find(cellIndex == 1);
    
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

    fileRE = fopen(['convertedClCdData/' profileNames{i} '_multiRe_Re.dat'],'w');
    fileCL = fopen(['convertedClCdData/' profileNames{i} '_multiRe_cl.dat'],'w');
    fileCD = fopen(['convertedClCdData/' profileNames{i} '_multiRe_cd.dat'],'w');

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
