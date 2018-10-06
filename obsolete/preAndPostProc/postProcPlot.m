% plot results obtained from postProc.m script

% input: postProcOutput directory
% output: on screen

clear all; close all; clc;
outputDir = "postProcOutput.12/";

fID_timeStr=fopen([outputDir "time.dat"]);
timeStr=textscan(fID_timeStr,'%s');
timeStr=timeStr{:};
nTime = length(timeStr);
time = zeros(nTime,1);
for i = 1:nTime
	time(i) = str2num(timeStr{i});
end
fclose(fID_timeStr);

turbineIDtoPlot = 1;
turbineDir = [outputDir "turbine" num2str(turbineIDtoPlot) "/"]
bladeDirList=dir([turbineDir "blade*"]); % get folder names
bladeDirListFlags = [bladeDirList.isdir]; % Get a logical vector that tells which is a directory.
bladeDirList=bladeDirList(bladeDirListFlags); % get those only directories
nBlades = length(bladeDirList);
plotLegendList = cell(nBlades,1);

for k = 1 : nBlades
	blade(k).name = bladeDirList(k).name;
	blade(k).nameTR = ["kanat" num2str(k)];
	plotLegendList{k} = blade(k).name;
	plotLegendListTR{k} = blade(k).nameTR;
	blade(k).ID = str2num(blade(k).name(6:end));
end

turbineResults = importdata([turbineDir "turbine" num2str(turbineIDtoPlot) ".dat"],' ',1).data;
turbineCp = turbineResults(:,1);
turbineCt = turbineResults(:,2);
bladeCp = turbineResults(:,3:2+nBlades);
bladeCt = turbineResults(:,3+nBlades:2+2*nBlades);

qbTurbineResults = importdata([turbineDir "qbTurbine" num2str(turbineIDtoPlot) ".dat"],' ',1).data;
qbTurbineCp = qbTurbineResults(1);
qbTurbineCt = qbTurbineResults(2);

qbBladeResults = importdata([turbineDir "qbTurbine" num2str(turbineIDtoPlot) ".Blade.dat"],' ',1).data;
qbBlade.radialPos   = qbBladeResults(:,1);
qbBlade.Fa          = qbBladeResults(:,2);
qbBlade.Fn          = qbBladeResults(:,3);
qbBlade.inflowAngle = qbBladeResults(:,4);
qbBlade.aoa         = qbBladeResults(:,5);
qbBlade.twist       = qbBladeResults(:,6);
qbBlade.re          = qbBladeResults(:,7);
qbBlade.cl          = qbBladeResults(:,8);
qbBlade.cd          = qbBladeResults(:,9);
qbBlade.vLoc        = qbBladeResults(:,10);
qbBlade.cn          = qbBladeResults(:,11);
qbBlade.ct          = qbBladeResults(:,12);

turbineConstants = importdata([turbineDir "constantsTurbine" num2str(turbineIDtoPlot) ".dat"],' ',1).data;
rotorRadius =  turbineConstants(1);
rotationDir =  turbineConstants(2);
velInf =       turbineConstants(3);
TSR =          turbineConstants(4);

% to plot blade1's desired positions
positionOfInterest = 180; % degrees
rotorSpeedDegPerSec = (velInf * TSR / rotorRadius)*180/pi;
rotorSecPerRev = 360 * 1 / rotorSpeedDegPerSec;
nRevolutions = 1+ floor ( ( max(time)-min(time) ) / rotorSecPerRev );
timeStamps = zeros(nRevolutions,1);
timeStamps(1) = min(time) + positionOfInterest / 360 * rotorSecPerRev;
for i=2:length(timeStamps)
	timeStamps(i) = timeStamps(i-1) + rotorSecPerRev;
end

clear qbBladeResults qbTurbineResults turbineConstants;

% --------------- TIME STEPPING RESULTS ---------------

bladePlotStyleList = {"b-o" "b-x" "b-s"};
bladePlotStyleList2 = {"b-" "b--" "b-."};
bladePlotStyleList3 = {"k-o" "k-x" "k-s"};
bladePlotStyleList4 = {"k-" "r-" "g-"};

plotColorsList = {[0 0 1] [0 0.5 0.5] [0.25 0.25 0.5]};

th = linspace(0,2*pi,100); %% to draw circle

plotForceDistr   = 0;
plotParameters   = 1;
plotCoeffs       = 1;
plot2Drotor      = 0; % should be 1 alone for real time rotation
plotPerformance  = 0;

for i =  nTime-80:nTime % [nTime-10 nTime-5 nTime]%
	figNum = 0;
	instant = timeStr{i};
	
	for b = 1:nBlades
		currentBladeDir = [turbineDir blade(b).name "/"];
		
		blade(b).data = importdata([currentBladeDir "tfDataOnElmCentrs." instant ".dat"],' ',1).data;
		blade(b).radialPosition = blade(b).data(:,1);
		blade(b).x              = blade(b).data(:,2);
		blade(b).y              = blade(b).data(:,3);
		blade(b).z              = blade(b).data(:,4);
		blade(b).aoa            = blade(b).data(:,5);
		blade(b).aoaGeom        = blade(b).data(:,6);
		blade(b).reynolds       = blade(b).data(:,7);
		blade(b).relVel         = blade(b).data(:,8);
		blade(b).fnMagPerSpan   = blade(b).data(:,9);
		blade(b).fxPerSpan      = blade(b).data(:,10);
		blade(b).fnVectorX      = blade(b).data(:,11);
		blade(b).fnVectorY      = blade(b).data(:,12);
		blade(b).liftCoeff      = blade(b).data(:,13);
		blade(b).dragCoeff      = blade(b).data(:,14);
		
		bladeInterp(b).data = importdata([currentBladeDir "tfDataOnInterpPts." instant ".dat"],' ',1).data;
		bladeInterp(b).radialPosition = bladeInterp(b).data(2:end-1,1);
		bladeInterp(b).errFn          = bladeInterp(b).data(2:end-1,2);
		bladeInterp(b).errFa          = bladeInterp(b).data(2:end-1,3);
		
	end
	%	minnak guş gaffası little birdy head

	if plotForceDistr == 1
		
		plotWindowSize = [0 0 10 10];
		
		%% PLOT NORMAL FORCE PER SPAN
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.Fn,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).fnMagPerSpan,bladePlotStyleList4{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","southeast");
		%title(["Normal kuvvet/uzunluk., time = " instant]);
		xlabel("Radyal pozisyon [m]");
		ylabel("Normal kuvvet/uzunluk [N/m]");
		title(["t = " num2str(time(i))]);
		ylim([600 1000])
		
		print(hFig,[outputDir "fnDist_" num2str(i) ".png"],"-dpng","-r0")
		
		%% PLOT AXIAL FORCE PER SPAN
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.Fa,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).fxPerSpan,bladePlotStyleList4{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","southeast");
		%title(["Axial Force per Length  Distr., time = " instant]);
		xlabel("Radyal pozisyon [m]");
		ylabel("Eksenel kuvvet/uzunluk [N/m]");
		title(["t = " num2str(time(i))]);
		
		%print(hFig,[outputDir "faDist_" num2str(i) ".png"],"-dpng","-r0")
		
		%%% PLOT NORMAL FORCE ERROR
		
		%subplotNum=subplotNum+1;
		%subplot(figureSize(1),figureSize(2),subplotNum);
		%for b = 1:nBlades
			%bar(bladeInterp(1).radialPosition,bladeInterp(1).errFn,...
			%'FaceColor',plotColorsList{b},'EdgeColor',"white");
			%hold on;
		%end
		%hold off;
		%grid on;
		%legend(plotLegendListTR{:});
		%title(["Normal force error., time = " instant]);
		%xlabel("Radial Position [m]");
		%ylabel("Relative error");
		%ylim([-1 1]);
		
		%%% PLOT AXIAL FORCE ERROR
		
		%subplotNum=subplotNum+1;
		%subplot(figureSize(1),figureSize(2),subplotNum);
		%for b = 1:nBlades
			%bar(bladeInterp(1).radialPosition,bladeInterp(1).errFa,...
			%'FaceColor',plotColorsList{b},'EdgeColor',"white");
			%hold on;
		%end
		%hold off;
		%grid on;
		%legend(plotLegendListTR{:});
		%title(["Axial force error., time = " instant]);
		%xlabel("Radial Position [m]");
		%ylabel("Relative error");
		%ylim([-1 1]);
		
	end
	
	if plotParameters == 1
		
		plotWindowSize = [0 0 10 10];
		
		%% PLOT LOCAL VELOCITY DISTRIBUTON
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.vLoc,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).relVel,bladePlotStyleList4{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","southeast");
		title(["t = " num2str(time(i))]);
		xlabel("Radyal Pozisyon [m]");
		ylabel("Bagil hiz [m/s]");
		ylim([20 80])
		
		%print(hFig,[outputDir "vRel_" num2str(i) ".png"],"-dpng","-r0")
		
		%% PLOT REYNOLDS DISTRIBUTON
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.re,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).reynolds,bladePlotStyleList{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","southeast");
		%title(["Reynolds number, time = " instant]);
		xlabel("Radyal Pozisyon [m]");
		ylabel("Reynolds Sayisi");
		
		print(hFig,"case2_graphs_re.png","-dpng","-r0")
		
		%% PLOT ANGLE OF ATTACK DISTRIBUTON
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.aoa,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).aoa,bladePlotStyleList4{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","northeast");
		title(["t = " num2str(time(i))]);
		xlabel("Radyal Pozisyon [m]");
		ylabel("Hucum acisi [deg]");
		
		ylim([4 8])
		
		print(hFig,[outputDir "aoa_" num2str(i) ".png"],"-dpng","-r0")
		
	end
	
	if plotCoeffs == 1
		
		plotWindowSize = [0 0 10 10];
		
		%% PLOT LIFT COEFF DISTRIBUTON
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.cl,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).liftCoeff,bladePlotStyleList4{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","northeast");
		title(["t = " num2str(time(i))]);
		xlabel("Radyal Pozisyon [m]");
		ylabel("Kaldirma Katsayisi, c_l");
		ylim([1 1.5])
		
		print(hFig,[outputDir "cl_" num2str(i) ".png"],"-dpng","-r0")
		
		%% PLOT DRAG COEFF DISTRIBUTON
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		plot(qbBlade.radialPos,qbBlade.cd,'b.');
		hold on;
		for b = 1:nBlades
			plot(blade(b).radialPosition,blade(b).dragCoeff,bladePlotStyleList{b});
		end
		hold off;
		grid on;
		legend("QBlade",plotLegendListTR{:},"location","northeast");
		%title(["Drag coefficient, time = " instant]);
		xlabel("Radyal Pozisyon [m]");
		ylabel("Direnc Katsayisi, C_d");
		
		print(hFig,"case2_graphs_cd.png","-dpng","-r0")
		
	end
	
	if plot2Drotor == 1
		
		plotWindowSize = [0 0 10 10];
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		if (exist("fig_blade")|exist("fig_vectors"))
			for b = 1:nBlades
				delete(fig_blade{b});
				delete(fig_vectors{b});
			end
			delete(plottedCirlce);
		end
		grid on;
		hold on;
		for b = 1:nBlades
			fig_blade{b} = plot(blade(b).y,blade(b).z,bladePlotStyleList4{b});
			fig_vectors{b} = quiver(blade(b).y,blade(b).z,...
			blade(b).fnVectorX ,blade(b).fnVectorY,'HandleVisibility','off');
			%quiver(0,0,unitX*rotorRadius,unitY*rotorRadius,'k','AutoScale','off');
			%quiver(0,0,unitNormalX*rotorRadius/2,unitNormalY*rotorRadius/2,'c','AutoScale','on');
		end
		plottedCirlce = plot(rotorRadius*cos(th),rotorRadius*sin(th),"k--");
		hold off;
		title(["t = " num2str(time(i))]);
		legend("Kanat 1","Kanat 2","Kanat 3");
		%xlabel("Alanin y-ekseni [m]");
		%ylabel("Alanin z-ekseni [m]");
		xlim([-rotorRadius*1.1 rotorRadius*1.1])
		ylim([-rotorRadius*1.1 rotorRadius*1.1])
		axis("square");
		%if i == nTime
			%dt = 0;
		%else
			%dt=time(i+1)-time(i);
		%end
		
		print(hFig,[outputDir "rotor2d_" num2str(i) ".png"],"-dpng","-r0")
	end

	pause(0)
end % end of time loop

% --------------- INSTANT RESULTS --------------- 

if plotPerformance == 1;
	
	plotWindowSize = [0 0 10 10];
	
	% PLOT TURBINE CP
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
	plot([0 time(nTime)],[16/27 16/27],'b-.');
	hold on;
	plot([0 time(nTime)],[qbTurbineCp qbTurbineCp],'b--');
	plot(time,turbineCp,'b-');
	hold off;
	grid on;
	%ylim([0 1]);
	legend("Betz limiti","QBlade","Hesaplanan","location","northeast");
	%title(["Turbine Power coefficient, time = " instant]);
	xlabel("Zaman [s]");
	ylabel("Turbin Guc Katsayisi, c_p");
		
		print(hFig,"case2_graphs_turbineCp.png","-dpng","-r0")
		
	% PLOT BLADE CPs
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
	hold on;
	for b = 1:nBlades
		plot(time,bladeCp(:,b),bladePlotStyleList2{b});
	end
	hold off;
	grid on;
	xlim([60 100]);
	legend(plotLegendListTR{:},"location","northeast");
	%title(["Blade Thrust coefficients, time = " instant]);
	xlabel("Zaman [s]");
	ylabel("Kanat Guc Katsayisi, c_p");
		
		print(hFig,"case2_graphs_bladeCp.png","-dpng","-r0")
		
	% PLOT TURBINE CT
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
	plot([0 time(nTime)],[qbTurbineCt qbTurbineCt],'b--');
	hold on;
	plot(time,turbineCt,'b-');
	hold off;
	grid on;
	%ylim([0 1]);
	legend("QBlade","Hesaplanan","location","northeast");
	%title(["Turbine Thrust coefficient, time = " instant]);
	xlabel("Zaman [s]");
	ylabel("Turbin Itki Katsayisi, c_t");
		
		print(hFig,"case2_graphs_turbineCt.png","-dpng","-r0")
		
	% PLOT BLADE CTs
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
	hold on;
	for b = 1:nBlades
		plot(time,bladeCt(:,b),bladePlotStyleList2{b});
	end
	hold off;
	grid on;
	xlim([60 100]);
	legend(plotLegendListTR{:},"location","northeast");
	%title(["Blade Thrust coefficients, time = " instant]);
	xlabel("Zaman [s]");
	ylabel("Kanat Itki Katsayisi, c_t");
		
		print(hFig,"case2_graphs_bladeCt.png","-dpng","-r0")
		
		nLastTimeStamps = 5;
		startTime = timeStamps(end-nLastTimeStamps);
		b=1; % which blade
		[minDifference startTimeIndex] = min(abs(time-startTime));
		
				close all
		
		% PLOT BLADE 1 CP WITH TIMESTAMPS
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		hold on;
		plot(time(startTimeIndex:nTime),bladeCp(startTimeIndex:nTime,b),bladePlotStyleList2{b});
		for ts = length(timeStamps)-nLastTimeStamps:length(timeStamps) 
			plot([timeStamps(ts) timeStamps(ts)],[min(bladeCp(startTimeIndex:nTime,b))/1.1 max(bladeCp(startTimeIndex:nTime,b))*1.1],"k--");
		end
		hold off;
		grid on;
		xlim([startTime time(end)])
		legend(plotLegendListTR{b},"Saat 6 pozisyonu","location","northeast");
		%title(["Blade Power coefficients, time = " instant]);
		xlabel("Zaman [s]");
		ylabel("Kanat Guc Katsayisi, c_p");
		
		print(hFig,"case2_graphs_blade1Cp.png","-dpng","-r0")
		
		% PLOT BLADE 1 CT WITH TIMESTAMPS
		
		figNum = figNum + 1;
		hFig = figure(figNum);
		set(hFig,'PaperUnits','centimeters','PaperPosition',plotWindowSize)
		
		hold on;
		plot(time(startTimeIndex:nTime),bladeCt(startTimeIndex:nTime,b),bladePlotStyleList2{b});
		for ts = length(timeStamps)-nLastTimeStamps:length(timeStamps) 
			plot([timeStamps(ts) timeStamps(ts)],[min(bladeCt(startTimeIndex:nTime,b))/1.1 max(bladeCt(startTimeIndex:nTime,b))*1.1],"k--");
		end
		hold off;
		grid on;
		xlim([startTime time(end)])
		legend(plotLegendListTR{b},"Saat 6 pozisyonu","location","northeast");
		%title(["Blade Power coefficients, time = " instant]);
		xlabel("Zaman [s]");
		ylabel("Kanat Itki Katsayisi, c_t");
		
		print(hFig,"case2_graphs_blade1Ct.png","-dpng","-r0")
		

end
	
	
	%% plot on blades visually (3D)
	%figure(2);
	%plot3(x,y,z,'-o');
	%hold on;
	%for k=1:nElements
		%t = linspace(0,2*pi,100);
		%zCirc = radialPosition(k)*cos(t);
		%yCirc = radialPosition(k)*sin(t);
		%xCirc = 0*t;
		%plot3(xCirc,yCirc,zCirc,'k-');
	%end
	%%quiver3(x,y,z,fxPerSpan,fyPerSpan,fzPerSpan);
	%hold off;
	%xlim([-rotorRadius rotorRadius])
	%ylim([-rotorRadius rotorRadius])
	%zlim([-rotorRadius rotorRadius])
	%title(["time = " num2str(instant)]);

%function orderPlot(figNum,subplotNum)
	%maxFigureSize = [3,3]; 
	%subplotNum = subplotNum + 1;
	%if (figNum == 0 | subplotNum > maxFigureSize(1)*maxFigureSize(2))
		%figNum = figNum + 1;
		%figure(figNum, 'position',[1 1 1300 700]);
		%subplotNum = 1;
	%end
	%subplot(maxFigureSize(1),maxFigureSize(2),subplotNum);
%end
