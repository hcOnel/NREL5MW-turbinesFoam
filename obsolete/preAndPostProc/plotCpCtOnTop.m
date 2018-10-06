clear all; close all; clc;

outputDirs = {"postProcOutput.10_2/" "postProcOutput.12/"};
plotStyleList = {"b-." "b-"};

figure(1,'position',[1 1 700 700]);
subplot(2,2,1);
hold on;
for i=1:length(outputDirs)
	fID_timeStr=fopen([outputDirs{i} "time.dat"]);
	timeStr=textscan(fID_timeStr,'%s');
	timeStr=timeStr{:};
	nTime = length(timeStr);
	time = zeros(nTime,1);
	for j = 1:nTime
		time(j) = str2num(timeStr{j});
	end
	fclose(fID_timeStr);
	
	turbineIDtoPlot = 1;
	turbineDir = [outputDirs{i} "turbine" num2str(turbineIDtoPlot) "/"]
	
	turbineResults = importdata([turbineDir "turbine" num2str(turbineIDtoPlot) ".dat"],' ',1).data;
	turbineCp = turbineResults(:,1);
	turbineCt = turbineResults(:,2);
	plot(time,turbineCp,plotStyleList{i})
	xlabel("Zaman, s");
	ylabel("Turbin Guc Katsayisi, c_p");
end
grid on

legend("Durum 1","Durum 2")
