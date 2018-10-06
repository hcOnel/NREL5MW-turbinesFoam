% compare turbinesFoam and qblade results

% input: postProcessing folder obtained by turbinesFoam and analysis results of qblade
% output: postProcOutput directory
% v02: multiple blades of a single turbine

clear all; close all; clc;
rho = 1.225; % air density, kg/m^3

turbinesFoamCaseDir = "postProcessing.08/";
turbinesSubdir	=			 [turbinesFoamCaseDir "turbines/0/"];
bladesSubdir	=		[turbinesFoamCaseDir "actuatorLines/0/"];
elementsSubdir	=[turbinesFoamCaseDir "actuatorLineElements/0/"];
qbSubdir = "rotorBEMsimulation_noTipLoss/";

%% IMPORT QBLADE RESULTS -----------------------------

qbTSRvalues = 1:0.5:10;

qbAxialForceData	= importdata([qbSubdir "Normal force Fn vs r.txt"],' ',3);
qbNormalForceData	= importdata([qbSubdir "Tangential force Ft vs r.txt"],' ',3);
qbInflowAngleData	= importdata([qbSubdir "Inflow Angle Phi vs r.txt"],' ',3);
qbAoAData			= importdata([qbSubdir "Angle of attack alpha vs r.txt"],' ',3);
qbTwistData			= importdata([qbSubdir "Blade twist angle theta vs r.txt"],' ',3);
qbReData			= importdata([qbSubdir "Reynolds number vs r.txt"],' ',3);
qbClData			= importdata([qbSubdir "Lift coeff vs r.txt"],' ',3);
qbCdData			= importdata([qbSubdir "Drag coeff vs r.txt"],' ',3);
qbVlocData			= importdata([qbSubdir "Local inflow speed vs r.txt"],' ',3);
qbCnData			= importdata([qbSubdir "Axial blade force coeff cn vs r.txt"],' ',3);
qbCtData			= importdata([qbSubdir "Tangential blade force coeff ct vs r.txt"],' ',3);
qbTurbineCpData		= importdata([qbSubdir "Power coeff vs TSR.txt"],' ',3);
qbTurbineCtData		= importdata([qbSubdir "Thrust coeff vs TSR.txt"],' ',3);

%% TURBINES -----------------------------

turbinesList=dir([turbinesSubdir "turbine*.csv"]); % get file names
turbinesList=struct2cell(turbinesList);
turbinesList=turbinesList(1,:)'; % get file names only as strings
nTurbines = length(turbinesList) % get number of files

time = (importdata([turbinesSubdir turbinesList{1}],',',1)).data(:,1);
nTime = length(time);
timeStr = cell(nTime,1);
for i = 1:nTime
	timeStr{i} = num2str(time(i),"%010.5f");
end

rotorRadiusData = [63]; % value in the fvOptions file
rotationDirectionData = [1]; % CW: 1, CCW: -1 (faced from upwind side, according to X direction)
velInfData = [11.4]; % free stream velocity, m/s
TSRdata = [7];
bladeAzimuthalOffset = [0 120 240]; % in degrees, !!! relative to the rotation axis! so write whatever is written on fvOptions

%IMPORTANT: blade visualization is seen from the downwind side!

outputDir = "postProcOutput.08/";
if isdir(outputDir)
	unix(["rm -r " outputDir]);
end
mkdir(outputDir);

fID_time = fopen([outputDir "time.dat"],"w");
fprintf(fID_time,"%s \n",timeStr{:});
fclose(fID_time);

for t=1:nTurbines
	locTurbines_ = strfind(turbinesList{t},'.'); % get locations of deliminators
	turbineNumber = turbinesList{t}(8:locTurbines_-1);
	turbineID = str2num(turbineNumber); % since it starts from 1
	
	currentTurbineDir = [outputDir "turbine" turbineNumber "/"];
	mkdir(currentTurbineDir);
	
	fID_constantsTurbine = fopen([currentTurbineDir "constantsTurbine" turbineNumber ".dat"],"w");
	fID_turbineData = fopen([currentTurbineDir "turbine" turbineNumber ".dat"],"w");
	fID_qbTurbineData = fopen([currentTurbineDir "qbTurbine" turbineNumber ".dat"],"w");
	fID_qbBladeData = fopen([currentTurbineDir "qbTurbine" turbineNumber ".Blade.dat"],"w");
	
	fprintf(fID_constantsTurbine,[repmat("%15s ",[1,4]) "\n"],"[1]rotorRadius","[2]rotationDir","[3]velInf","[4]TSR")
	fprintf(fID_qbTurbineData,[repmat("%15s ",[1,2]) "\n"],"[1]qbCp","[2]qbCt")
	fprintf(fID_qbBladeData,[repmat("%15s ",[1,12]) "\n"],"[1]radialPos","[2]Fa","[3]Fn",...
		"[4]inflowAngle","[5]aoa","[6]twist","[7]re","[8]cl","[9]cd","[10]vLoc","[11]cn(ca)","[12]ct(cn)")
	
	turbineData = (importdata([turbinesSubdir turbinesList{t}],',',1)).data; % has rows as much as nTime
	rotationDirection = rotationDirectionData(turbineID);
	velInf = velInfData(turbineID);
	TSR = TSRdata(turbineID);
	rotorRadius = rotorRadiusData(turbineID);
	
	fprintf(fID_constantsTurbine,"%15e %15i %15e %15e \n",rotorRadius,rotationDirection,velInf,TSR)
	
	rotorCircumference = 2*pi*rotorRadius; % m
	tipSpeed = TSR*velInf; % m/s
	rotorSpeed = 2*pi*tipSpeed/rotorCircumference; % rad/s
	
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
	
	%% GET QBLADE RESULTS FOR TURBINE TSR -----------------------------
	
	TSRindex = find(qbTSRvalues==TSR)*2-1;
	TSRindexTurbine = find(qbTSRvalues==TSR);
	
	qbRadialPos		= qbAoAData			.data(:,TSRindex); % get the list of radial positions from AoA, assuming it remains the same for all parameters
	
	qbAxialForce	= qbAxialForceData	.data(:,TSRindex:TSRindex+1);
	qbNormalForce	= qbNormalForceData	.data(:,TSRindex:TSRindex+1);
	qbInflowAngle	= qbInflowAngleData	.data(:,TSRindex:TSRindex+1);
	qbAoA			= qbAoAData			.data(:,TSRindex:TSRindex+1);
	qbTwist			= qbTwistData		.data(:,TSRindex:TSRindex+1);
	qbRe			= qbReData			.data(:,TSRindex:TSRindex+1);
	qbCl			= qbClData			.data(:,TSRindex:TSRindex+1);
	qbCd			= qbCdData			.data(:,TSRindex:TSRindex+1);
	qbVloc			= qbVlocData		.data(:,TSRindex:TSRindex+1);
	qbCn			= qbCnData			.data(:,TSRindex:TSRindex+1);
	qbCt			= qbCtData			.data(:,TSRindex:TSRindex+1);
	qbTurbineCp		= qbTurbineCpData	.data(TSRindexTurbine,2);
	qbTurbineCt		= qbTurbineCtData	.data(TSRindexTurbine,2);
	
	fprintf(fID_qbBladeData,[repmat("%15e ",[1,12]) "\n"],[qbRadialPos';...
	qbAxialForce(:,2)';qbNormalForce(:,2)';qbInflowAngle(:,2)';qbAoA(:,2)';...
	qbTwist(:,2)';qbRe(:,2)';qbCl(:,2)';qbCd(:,2)';qbVloc(:,2)';qbCn(:,2)';qbCt(:,2)'])
	fprintf(fID_qbTurbineData,[repmat("%15e ",[1,2]) "\n"],qbTurbineCp,qbTurbineCt)
	
	rotorAngleDeg = turbineData(:,2);
	rotorAngle = turbineData(:,2)*pi/180*rotationDirection;
	cpTurbine = turbineData(:,4);
	ctTurbine = turbineData(:,5);
	
	clear turbineData
	
	%% BLADES -----------------------------
	
	bladesList=dir([bladesSubdir "turbine" mat2str(turbineID) "*.csv"]); % get file names
	bladesList=struct2cell(bladesList);
	bladesList=bladesList(1,:)'; % get file names only as strings
	nBlades = length(bladesList) % get number of files
	
	manualTurbineCp = zeros(nTime,1);
	manualTurbineCt = zeros(nTime,1);
	bladeAngularPosDeg = zeros(nTime,nBlades);
	manualBladeCp = zeros(nTime,nBlades);
	manualBladeCt = zeros(nTime,nBlades);
	
	for b=1:nBlades
		locBlades_ = strfind(bladesList{b},'.'); % get locations of deliminators
		bladeNumber = bladesList{b}(locBlades_(1)+6:locBlades_(2)-1)
		bladeID = str2num(bladeNumber); % since it starts from 1
		
		currentBladeDir = [currentTurbineDir "blade" bladeNumber "/"];
		mkdir(currentBladeDir);
		
		bladeData = (importdata([turbinesSubdir turbinesList{t}],',',1)).data;
		
		%% ELEMENTS -----------------------------
		
		elementsList=dir([elementsSubdir "turbine" mat2str(turbineID) ".blade" mat2str(bladeID) "*.csv"]); % get file names
		elementsList=struct2cell(elementsList);
		elementsList=elementsList(1,:)'; % get file names only as strings
		nElements = length(elementsList) % get number of files
		bladeElementData = cell(nTime,nElements); % initialize cell that will store everything
		
		for e=1:nElements
			loc_ = strfind(elementsList{e},'.'); % get locations of deliminators
			elementNumber = elementsList{e}(loc_(2)+8:loc_(3)-1);
			elementID = str2num(elementNumber)+1; % since it starts from 0
			elementData = (importdata([elementsSubdir elementsList{e}],',',1)).data;
			for j=1:nTime
				bladeElementData{j,elementID}=elementData(j,:);
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

		%% CALCULATIONS -------------------------------

		radialPosition = zeros(nElements,1);
		x = zeros(nElements,1);
		y = zeros(nElements,1);
		z = zeros(nElements,1);
		relVel = zeros(nElements,1);
		reynolds = zeros(nElements,1);
		aoa = zeros(nElements,1);
		aoaGeom = zeros(nElements,1);
		liftCoeff = zeros(nElements,1);
		dragCoeff = zeros(nElements,1);
		fnMagPerSpan = zeros(nElements,1);
		%fnMagPerSpan2 = zeros(nElements,1); % to compare
		fnVectorX = zeros(nElements,1);
		fnVectorY = zeros(nElements,1);
		endEffectFactor = zeros(nElements,1);
		%elementCp = zeros(nElements,1);
		%elementCt = zeros(nElements,1);

		% --- For manual calculations
		windPower = 0.5*rho*(pi*rotorRadius^2)*velInf^3;
		windForce = 0.5*rho*(pi*rotorRadius^2)*velInf^2;
		fvOptionsSectionData = load("turbine1.blade1.tfSectionGeometry");
		chordElement = fvOptionsSectionData(:,4);
		pitchElement = fvOptionsSectionData(:,6);
		spanElement = zeros(nElements,1);
		for e=1:nElements
			spanElement(e)=fvOptionsSectionData(e+1,2)-fvOptionsSectionData(e,2);
		end
		linearSpanPoints = linspace(0,rotorRadius); % for interpolation in use of error
		%fLiftElement = zeros(nElements,1);
		%fDragElement = zeros(nElements,1);
		%fAxialElement = zeros(nElements,1);
		%fNormalElement = zeros(nElements,1);

										bladePlotSwitch = 1;
										turbine2DvisualPlotSwitch = 0;
										turbinePlotSwitch = 0;

		for i=1:nTime
			instant=time(i);
			
			fID_tfDataOnElmCentrs = fopen([currentBladeDir "tfDataOnElmCentrs." timeStr{i} ".dat"],"w");
			fID_tfDataOnInterpPts = fopen([currentBladeDir "tfDataOnInterpPts." timeStr{i} ".dat"],"w");
			
			fprintf(fID_tfDataOnElmCentrs,[repmat("%15s ",[1,14]) "\n"],"[1]radialPos","[2]xCoord","[3]yCoord","[4]zCoord",...
				"[5]aoa","[6]aoaGeom","[7]reynolds","[8]relVelMag","[9]fnMagPerSpn","[10]faMagPerSpn",...
				"[11]fnVectorX","[12]fnVectorY","[13]cl","[14]cd")
			fprintf(fID_tfDataOnInterpPts,[repmat("%15s ",[1,3]) "\n"],"[1]radialPos","[2]errFn","[3]errFa")
				
			% create the unit vector of blade direction, add 90 degrees since 
			% the positive vertical is assumed zero degrees in turbinesFoam
			% also use rotation direction
			% v02: also add blade initial angle
			bladeInitialAngle = rotationDirection * bladeAzimuthalOffset(turbineID,bladeID)*pi/180;
			instantBladeAngle = bladeInitialAngle + (rotorAngle(i)+pi/2);
			instantNormalAngle = (instantBladeAngle+rotationDirection*pi/2);
			[unitX unitY] = pol2cart(instantBladeAngle,1);
			[unitNormalX unitNormalY] = pol2cart(instantNormalAngle,1);
			bladeUnitVector = [unitX unitY];
			bladeNormalUnitVector = [unitNormalX unitNormalY];
			
			manualBladeTorque = 0;
			manualBladeThrust = 0;
			
			for j=1:nElements
				
				% check table for index of what to plot before!
				radialPosition(j) = bladeElementData{i,j}(2)*rotorRadius;
				x(j) = bladeElementData{i,j}(3);
				y(j) = bladeElementData{i,j}(4);
				z(j) = bladeElementData{i,j}(5);
				relVel(j) = bladeElementData{i,j}(6);
				reynolds(j) = bladeElementData{i,j}(7);
				fx(j) = bladeElementData{i,j}(12);% *rho;
				fy(j) = bladeElementData{i,j}(13);% *rho;
				fz(j) = bladeElementData{i,j}(14);% *rho;
				fxPerSpan(j) = fx(j)/spanElement(j);
				fyPerSpan(j) = fy(j)/spanElement(j);
				fzPerSpan(j) = fz(j)/spanElement(j);
				liftCoeff(j) = bladeElementData{i,j}(10);
				dragCoeff(j) = bladeElementData{i,j}(11);
				fnMagPerSpan(j) = dot([fyPerSpan(j) fzPerSpan(j)],bladeNormalUnitVector); % component of the local force in the normal direction
				[fnVectorX(j) fnVectorY(j)] = pol2cart(instantNormalAngle,fnMagPerSpan(j));
				%fnMagPerSpan2(j)=norm([fyPerSpan(j) fzPerSpan(j)]);
				aoa(j) = bladeElementData{i,j}(8);
				aoaGeom(j) = bladeElementData{i,j}(9);
				relVel(j) = bladeElementData{i,j}(6);
				
					%%calculate forces manually
					%fLiftElement(j) = 0.5*rho*relVel(j)^2*liftCoeff(j)*chordElement(j)*spanElement(j);
					%fDragElement(j) = 0.5*rho*relVel(j)^2*dragCoeff(j)*chordElement(j)*spanElement(j);
					%axis2liftAngle =  (aoa(j)+pitchElement(j))*pi/180 ;
					%fAxialElement = fLiftElement*cos(axis2liftAngle) + fDragElement*sin(axis2liftAngle);
					%fNormalElement = fLiftElement*sin(axis2liftAngle) - fDragElement*cos(axis2liftAngle);
				
				manualBladeTorque = manualBladeTorque + fnMagPerSpan(j)*spanElement(j)*radialPosition(j);
				manualBladeThrust = manualBladeThrust + fxPerSpan(j)*spanElement(j);
				
			end % end of nElements
			
			manualBladeCp(i,b) = manualBladeTorque*rotorSpeed / windPower;
			manualBladeCt(i,b) = manualBladeThrust / windForce;
			bladeAngularPosDeg(i,b) = rotorAngleDeg(i) + bladeAzimuthalOffset(b);
			manualTurbineCp(i) = manualTurbineCp(i) + manualBladeCp(i,b);
			manualTurbineCt(i) = manualTurbineCt(i) + manualBladeCt(i,b);
			
			errFnFoam = interp1(radialPosition,fnMagPerSpan,linearSpanPoints,"linear","extrap");
			errFnQb = interp1(qbNormalForce(:,1),qbNormalForce(:,2),linearSpanPoints,"linear","extrap");
			errFn = (errFnFoam - errFnQb) ./ errFnQb;
			errFaFoam = interp1(radialPosition,fxPerSpan,linearSpanPoints,"linear","extrap");
			errFaQb = interp1(qbAxialForce(:,1),qbAxialForce(:,2),linearSpanPoints,"linear","extrap");
			errFa = (errFaFoam - errFaQb) ./ errFaQb;
			
			fprintf(fID_tfDataOnElmCentrs,[repmat("%15e ",[1,14]) "\n"],[radialPosition';x';y';z';aoa';aoaGeom';...
			reynolds';relVel';fnMagPerSpan';fxPerSpan;fnVectorX';fnVectorY';liftCoeff';dragCoeff']);
			fprintf(fID_tfDataOnInterpPts,[repmat("%15e ",[1,3]) "\n"],[linearSpanPoints;errFn;errFa]);

			fclose(fID_tfDataOnElmCentrs);
			fclose(fID_tfDataOnInterpPts);
			
		end % end of nTime
		
	end % end of nBlades
	
	fprintf(fID_turbineData,[repmat("%15s ",[1,3]) "\n"],"[1]cp","[2]ct","| cp(all blades) | ct(all blades) | in order")
	fprintf(fID_turbineData,[repmat("%15e ",[1,2+2*nBlades]) "\n"],[manualTurbineCp';manualTurbineCt';manualBladeCp';manualBladeCt'])
	
	fclose(fID_constantsTurbine);
	fclose(fID_turbineData);
	fclose(fID_qbTurbineData);
	fclose(fID_qbBladeData);
end % end of nTurbines
fclose("all");
disp("Done.")
