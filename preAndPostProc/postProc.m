% compare turbinesFoam and qblade results

% input: postProcessing folder obtained by turbinesFoam and analysis results of qblade
% output: on screen

clear all; close all; clc;

% WORKS FOR ONLY 1 TURBINE

rotorRadius = 63; % value in the fvOptions file
rotationDirection = 1; % CW: 1, CCW: -1 (faced from upwind side, according to X direction)
	%IMPORTANT: blade visualization is seen from the downwind side!
TSR = 7;
rho = 1.225; % air density, kg/m^3
velInf = 11.4; % free stream velocity, m/s
nBlades = 3;

rotorCircumference = 2*pi*rotorRadius; % m
tipSpeed = TSR*velInf; % m/s
rotorSpeed = 2*pi*tipSpeed/rotorCircumference; % rad/s

turbinesFoamCaseDir = "postProcessing/";
turbinesSubdir= [turbinesFoamCaseDir "turbines/0.8/"];
bladesSubdir= [turbinesFoamCaseDir "actuatorLines/0.8/"];
elementsSubdir= [turbinesFoamCaseDir "actuatorLineElements/0.8/"];
qBladeSubdir = "rotorBEMsimulation_noTipLoss/";

%% IMPORT QBLADE RESULTS -----------------------------

qBladeTSRvalues = 1:0.5:10;
TSRindex = find(qBladeTSRvalues==TSR)*2-1;

qBladeAxialForceData	= importdata([qBladeSubdir "Normal force Fn vs r.txt"],' ',3);
qBladeNormalForceData	= importdata([qBladeSubdir "Tangential force Ft vs r.txt"],' ',3);
qBladeInflowAngleData	= importdata([qBladeSubdir "Inflow Angle Phi vs r.txt"],' ',3);
qBladeAoAData			= importdata([qBladeSubdir "Angle of attack alpha vs r.txt"],' ',3);
qBladeTwistData			= importdata([qBladeSubdir "Blade twist angle theta vs r.txt"],' ',3);
qBladeReData			= importdata([qBladeSubdir "Reynolds number vs r.txt"],' ',3);
qBladeClData			= importdata([qBladeSubdir "Lift coeff vs r.txt"],' ',3);
qBladeCdData			= importdata([qBladeSubdir "Drag coeff vs r.txt"],' ',3);
qBladeVlocData			= importdata([qBladeSubdir "Local inflow speed vs r.txt"],' ',3);
qBladeCnData			= importdata([qBladeSubdir "Axial blade force coeff cn vs r.txt"],' ',3);
qBladeCtData			= importdata([qBladeSubdir "Tangential blade force coeff ct vs r.txt"],' ',3);

qBladeAxialForce	= qBladeAxialForceData	.data(:,TSRindex:TSRindex+1);
qBladeNormalForce	= qBladeNormalForceData	.data(:,TSRindex:TSRindex+1);
qBladeInflowAngle	= qBladeInflowAngleData	.data(:,TSRindex:TSRindex+1);
qBladeAoA			= qBladeAoAData			.data(:,TSRindex:TSRindex+1);
qBladeTwist			= qBladeTwistData		.data(:,TSRindex:TSRindex+1);
qBladeRe			= qBladeReData			.data(:,TSRindex:TSRindex+1);
qBladeCl			= qBladeClData			.data(:,TSRindex:TSRindex+1);
qBladeCd			= qBladeCdData			.data(:,TSRindex:TSRindex+1);
qBladeVloc			= qBladeVlocData		.data(:,TSRindex:TSRindex+1);
qBladeCn			= qBladeCnData			.data(:,TSRindex:TSRindex+1);
qBladeCt			= qBladeCtData			.data(:,TSRindex:TSRindex+1);

clear...
qBladeAxialForceData	...
qBladeNormalForceData	...
qBladeInflowAngleData	...
qBladeAoAData			...
qBladeTwistData			...
qBladeReData			...
qBladeClData			...
qBladeCdData			...
qBladeVlocData			...
qBladeCnData			...
qBladeCtData			;

%% TURBINES -----------------------------

turbinesList=dir([turbinesSubdir "*.csv"]); % get file names
turbinesList=struct2cell(turbinesList);
turbinesList=turbinesList(1,:)'; % get file names only as strings
nTurbines = length(turbinesList) % get number of files
turbineDataTemp2 = importdata([turbinesSubdir turbinesList{1}],',',1);
turbineDataTemp = turbineDataTemp2.data;

%	[1,1] = time
%	[1,2] = angle_deg
%	[1,3] = tsr
%	[1,4] = cp
%	[1,5] = cd
%	[1,6] = ct
%	[1,7] = cd_blade1
%	[1,8] = ct_blade1
%	[1,9] = cd_blade2
%	[1,10] = ct_blade2
%	[1,11] = cd_blade3
%	[1,12] = ct_blade3

time=turbineDataTemp(:,1);
nTime = length(time);
rotorAngle = turbineDataTemp(:,2)*pi/180*rotationDirection;
cpTurbine = turbineDataTemp(:,4);
ctTurbine = turbineDataTemp(:,5);

clear turbineDataTemp2 turbineDataTemp

%% BLADES -----------------------------

	% nothing yet

%% ELEMENTS -----------------------------

elementsList=dir([elementsSubdir "*.csv"]); % get file names
elementsList=struct2cell(elementsList);
elementsList=elementsList(1,:)'; % get file names only as strings
nElements = length(elementsList) % get number of files
bladeElementData = cell(nTime,nElements); % initialize cell that will store everything

for i=1:nElements
	loc_ = strfind(elementsList{i},'.'); % get locations of deliminators
	elementNumber = elementsList{i}(loc_(2)+8:loc_(3)-1);
	elementNumber = str2num(elementNumber)+1; % since it starts from 0
	elementDataTemp2 = importdata([elementsSubdir elementsList{i}],',',1);
	elementDataTemp = elementDataTemp2.data;
	for j=1:nTime
		bladeElementData{j,elementNumber}=elementDataTemp(j,:);
	end
end

%	[1,1] = time
%	[1,2] = root_dist
%	[1,3] = x
%	[1,4] = y
%	[1,5] = z
%	[1,6] = rel_vel_mag
%	[1,7] = Re
%	[1,8] = alpha_deg
%	[1,9] = alpha_geom_deg
%	[1,10] = cl
%	[1,11] = cd
%	[1,12] = fx
%	[1,13] = fy
%	[1,14] = fz
%	[1,15] = end_effect_factor

%% CALCULATIONS AND PLOTS -------------------------------

radialPosition = zeros(1,nElements);
x = zeros(1,nElements);
y = zeros(1,nElements);
z = zeros(1,nElements);
relVel = zeros(1,nElements);
reynolds = zeros(1,nElements);
aoa = zeros(1,nElements);
aoaGeom = zeros(1,nElements);
liftCoeff = zeros(1,nElements);
dragCoeff = zeros(1,nElements);
fNormalMagnitude = zeros(1,nElements);
%fNormalMagnitude2 = zeros(1,nElements); % to compare
fNormalVectorX = zeros(1,nElements);
fNormalVectorY = zeros(1,nElements);
endEffectFactor = zeros(1,nElements);

% --- For manual calculations
windPower = 0.5*rho*(pi*rotorRadius^2)*velInf^3;
windForce = 0.5*rho*(pi*rotorRadius^2)*velInf^2;
manCalcCp = zeros(nTime,1);
manCalcCt = zeros(nTime,1);
fvOptionsSectionData = load('fvOptionsBladeGeomInput.dat');
chordElement = fvOptionsSectionData(:,4);
pitchElement = fvOptionsSectionData(:,6);
spanElement = zeros(1,nElements);
for i=1:nElements
	spanElement(i)=fvOptionsSectionData(i+1,2)-fvOptionsSectionData(i,2);
end
%fLiftElement = zeros(1,nElements);
%fDragElement = zeros(1,nElements);
%fAxialElement = zeros(1,nElements);
%fNormalElement = zeros(1,nElements);

bladePlotSwitch = 1;
turbinePlotSwitch = 0;

for i=1:nTime
	instant=time(i);
	% create the unit vector of blade direction, add 90 degrees since 
	% the positive vertical is assumed zero degrees in turbinesFoam
	% also use rotation direction
	instantBladeAngle = (rotorAngle(i)+pi/2);
	instantNormalAngle = (instantBladeAngle+rotationDirection*pi/2);
	[unitX unitY] = pol2cart(instantBladeAngle,1);
	[unitNormalX unitNormalY] = pol2cart(instantNormalAngle,1);
	bladeUnitVector = [unitX unitY];
	bladeNormalUnitVector = [unitNormalX unitNormalY];
	manCalcTorque = 0;
	manCalcThrust = 0;

	for j=1:nElements

		% check table for index of what to plot before!
		radialPosition(j) = bladeElementData{i,j}(2)*rotorRadius;
		x(j) = bladeElementData{i,j}(3);
		y(j) = bladeElementData{i,j}(4);
		z(j) = bladeElementData{i,j}(5);
		relVel(j) = bladeElementData{i,j}(6);
		reynolds(j) = bladeElementData{i,j}(7);
		fx(j) = bladeElementData{i,j}(12);
		fy(j) = bladeElementData{i,j}(13);
		fz(j) = bladeElementData{i,j}(14);
		fxPerSpan(j) = fx(j)/spanElement(j);
		fyPerSpan(j) = fy(j)/spanElement(j);
		fzPerSpan(j) = fz(j)/spanElement(j);
		liftCoeff(j) = bladeElementData{i,j}(10);
		dragCoeff(j) = bladeElementData{i,j}(11);
		fNormalMagnitude(j) = dot([fyPerSpan(j) fzPerSpan(j)],bladeNormalUnitVector); % component of the local force in the normal direction
		[fNormalVectorX(j) fNormalVectorY(j)] = pol2cart(instantNormalAngle,fNormalMagnitude(j));
		%fNormalMagnitude2(j)=norm([fyPerSpan(j) fzPerSpan(j)]);
		aoa(j) = bladeElementData{i,j}(8);
		aoaGeom(j) = bladeElementData{i,j}(9);
		relVel(j) = bladeElementData{i,j}(6);
		
			%%calculate forces manually
			%fLiftElement(j) = 0.5*rho*relVel(j)^2*liftCoeff(j)*chordElement(j)*spanElement(j);
			%fDragElement(j) = 0.5*rho*relVel(j)^2*dragCoeff(j)*chordElement(j)*spanElement(j);
			%axis2liftAngle =  (aoa(j)+pitchElement(j))*pi/180 ;
			%fAxialElement = fLiftElement*cos(axis2liftAngle) + fDragElement*sin(axis2liftAngle);
			%fNormalElement = fLiftElement*sin(axis2liftAngle) - fDragElement*cos(axis2liftAngle);
		
		manCalcTorque = manCalcTorque + fNormalMagnitude(j)*spanElement(j)*radialPosition(j);
		manCalcThrust = manCalcThrust + fxPerSpan(j)*spanElement(j);
	end
	
	manCalcPower = manCalcTorque*rotorSpeed*nBlades;
	manCalcCp(i) = manCalcPower / windPower;
	manCalcCt(i) = manCalcThrust*nBlades / windForce;
	
	if bladePlotSwitch == 1
		
		figure(1, 'position',[1 1 1200 700]);
		figureSize = [2,3]; figNum = 0;
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % normal force
		plot(radialPosition,fNormalMagnitude,"b-");
		hold on;
		% qBlade results
		plot(qBladeNormalForce(:,1),qBladeNormalForce(:,2),"b.");
		%	minnak guş gaffası little birdy head
		hold off;
		title(["Normal Force Distro per length, time = " num2str(instant)]);
		legend("Fn","QBlade Fn",'location','northwest');
		xlabel("Radial position [m]");
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % axial force
		plot(radialPosition,fxPerSpan,"r-");
		hold on;
		% qBlade results
		plot(qBladeAxialForce(:,1),qBladeAxialForce(:,2),"r.");
		hold off;
		title(["Axial Force Distro per length, time = " num2str(instant)]);
		legend("Fa","QBlade Fa",'location','northwest');
		xlabel("Radial position [m]");
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % angle of attack
		plot(radialPosition,aoa,"b-");
		hold on;
		plot(radialPosition,aoaGeom,"r-");
		% qBlade results
		plot(qBladeAoA(:,1),qBladeAoA(:,2),"b.");
		hold off;
		title(["Angle of Attack, time = " num2str(instant)]);
		legend("AoA","AoA Geom","QBlade Alpha")
		xlabel("Radial position [m]");
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % reynolds
		plot(radialPosition,reynolds,"b-");
		hold on;
		% qBlade results
		plot(qBladeRe(:,1),qBladeRe(:,2),"b.");
		hold off;
		title(["Reynolds Local, time = " num2str(instant)]);
		legend("Re","QBlade Re",'location','northwest')
		xlabel("Radial position [m]");
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % relative velocity
		plot(radialPosition,relVel,"b-");
		hold on;
		% qBlade results
		plot(qBladeVloc(:,1),qBladeVloc(:,2),"b.");
		hold off;
		title(["Local velocity, time = " num2str(instant)]);
		legend("Vloc","QBlade Vloc",'location','northwest');
		xlabel("Radial position [m]");
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % cl cd
		plot(radialPosition,liftCoeff,"b-");
		hold on;
		plot(radialPosition,dragCoeff,"r-");
		% qBlade results
		plot(qBladeCl(:,1),qBladeCl(:,2),"b.");
		plot(qBladeCd(:,1),qBladeCd(:,2),"r.");
		%plot(radialPosition,fNormalMagnitude2,"kx");
		hold off;
		title(["Cl Cd, time = " num2str(instant)]);
		legend("Cl","Cd","QBlade Cl","QBlade Cd")
		xlabel("Radial position [m]");
		
		%figure(10, 'position',[1 1 500 500]);
		%plot(y,z,'ro');
		%hold on;
		%%quiverScale=0.001;
		%%quiver(y,z,fyPerSpan,fzPerSpan);
		%quiver(y,z,fNormalVectorX ,fNormalVectorY );
		%quiver(0,0,unitX*rotorRadius,unitY*rotorRadius,'k','AutoScale','off');
		%quiver(0,0,unitNormalX*rotorRadius/2,unitNormalY*rotorRadius/2,'c','AutoScale','on');
		%hold off;
		%xlim([-rotorRadius rotorRadius])
		%ylim([-rotorRadius rotorRadius])
		%title(["2D blade visualization, time = " num2str(instant)]);
		%legend("Blade element centers","Local normal force");
		%xlabel("Radial position [m]");
		%grid on;
		 
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
	end
	
	if turbinePlotSwitch == 1

		figure(2, 'position',[1 1 1200 700]);
		figureSize = [1,2]; figNum = 0;
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % power curve
		plot(time,cpTurbine);
		hold on;
		% qBlade results
		plot(time,manCalcCp,'r');
		hold off;
		title(["Power coefficient Cp, time = " num2str(instant)]);
		legend("Cp","Manual Cp")
		xlabel("Time step [s]");
		
		figNum=figNum+1;
		subplot(figureSize(1),figureSize(2),figNum); % thrust curve
		plot(time,ctTurbine);
		hold on;
		% qBlade results
		plot(time,manCalcCt,'r');
		hold off;
		title(["Thrust coefficient Cp, time = " num2str(instant)]);
		legend("Ct","Manual Ct")
		xlabel("Time step [s]");

	end
	pause(0.0);
end


