% converts qblade's blade sections data to turbineFoam's fvOptions input format (elements)

% input: bladeGeom.txt and bladeProfiles.txt files
% output: fvOptionsBladeGeomInput.dat


clear all;close all;clc;
bladeGeomData=load("bladeGeom.txt");
profilesFileID=fopen('bladeProfiles.txt');
bladeProfileNames=textscan(profilesFileID,'%s');
nQbladeSections=length(bladeGeomData(:,1));

nElements = nQbladeSections;
nSections = nElements + 1;
elementData = zeros(nSections,6);

% input (QBlade's output) format:
%        Radial Position [m]          Chord Length [m]               Twist [deg]     Pitch Axis Offset [m]  Thread Axis in [% chord]              Airfoil Name            360 Polar Name
qbRadius = bladeGeomData(:,1);
qbChord = bladeGeomData(:,2);
qbTwist = bladeGeomData(:,3);
qbPitchAxisOffset = -bladeGeomData(:,4);
qbThreadAxisPercent = bladeGeomData(:,5);

clear bladeGeomData;

qbChordMountDist = (qbPitchAxisOffset+qbChord.*qbThreadAxisPercent);
qbChordMountPercent = qbChordMountDist./qbChord;

% output (fvOptions) format: 
% elementData:
% axialDistance, radius, azimuth, chord, chordMount, pitch

chordMountDist = zeros(1,nSections);

for i=1:nSections
	if i==1
		lefti = i;
		righti = i;
	elseif i==nSections
		lefti = i-1;
		righti = i-1;
	else
		lefti = i-1;
		righti = i;
	end
	% axial distance is always zero in QBlade !!!
	radius=(qbRadius(lefti) + qbRadius(righti))/2;
	elementData(i,2)=radius;
	% azimuth is always zero in QBlade !!!
	chord = interp1(qbRadius,qbChord,radius,"cubic");
	elementData(i,4)=chord;
	chordMountDistTemp = interp1(qbRadius,qbChordMountDist,radius,"cubic");
	chordMountDist(i) = chordMountDistTemp;
	elementData(i,5)=chordMountDistTemp/chord;
	pitch = interp1(qbRadius,qbTwist,radius,"cubic");
	elementData(i,6)=pitch;
end

fid = fopen("fvOptionsBladeGeomInput.dat",'w');

fprintf('nElements %i;\n', nElements);
%fprintf('( %e %e %e %e %e %e )\n', elementData.');
fprintf(fid,'%e %e %e %e %e %e \n', elementData.');
fclose(fid);

figure;
plot(elementData(:,2),elementData(:,4),'b-');
hold on;
plot(qbRadius,qbChord,'bx');

plot(elementData(:,2),chordMountDist,'r-');
plot(qbRadius,qbChordMountDist,'rx');

%plot(elementData(:,2),elementData(:,5),'k-');
%plot(qbRadius,qbChordMountPercent,'kx');

plot(elementData(:,2),elementData(:,6),'k-');
plot(qbRadius,qbTwist,'kx');

legend('chord','qbChord','chordMountDist','qbChordMount','pitch','qbTwist');




