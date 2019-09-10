%All units are SI, always.

CylinderHeight = 5e-2;
CylinderDiameter = 2.54e-2;
CylinderMass = 1;
CylinderAxialGridPoints = 5;
CylinderRadialGridPoints = 5;

RotorRadius2 = 10e-2;
RotorRadius3 = 15e-2;

RotorPosition = [2 0 0];

TMLength = 30e-2;
TMDiameter = 40e-2;
TMMass = 40;
TMAxialGridPoints = 20;
TMRadialGridPoints = 20;


%Make a rotor cylinder
Cylinder = genPointMassAnnlSheet(CylinderMass, 0, CylinderDiameter/2, ...
		CylinderHeight, ...
		CylinderAxialGridPoints, CylinderRadialGridPoints);
Cylinder = rotatePMArray(Cylinder, pi/2, [0 1 0]);	


%Make the rotor arrangement
R2 = translatePMArray(Cylinder, [RotorRadius2 0 0]);
R3 = translatePMArray(Cylinder, [RotorRadius3 0 0]);

Rotor = [
		R2;
		rotatePMArray(R2, pi, [0 0 1]);
		R3;
		rotatePMArray(R3, 2*pi/3, [0 0 1]);
		rotatePMArray(R3, 4*pi/3, [0 0 1]);
	];

%Make a LIGO TM
TM = genPointMassAnnlSheet(TMMass, 0, TMDiameter/2, TMLength,...
		TMAxialGridPoints, TMRadialGridPoints);

%Spin the rotor and compute the torques

AngleForceTorque = [];

for angle = 0:0.05:(2*pi)*1

	angle  %For humans

	RotatedRotor = rotatePMArray(Rotor, angle, [0 0 1]);
	TranslatedRotatedRotor = translatePMArray(RotatedRotor, RotorPosition);

	displayPoints(TM, TranslatedRotatedRotor); 
	[Force Torque ] = pointMatrixGravity(TM, TranslatedRotatedRotor);

	AngleForceTorque = [AngleForceTorque; angle Force Torque];
end

%%

figure(2)
l=plot(AngleForceTorque(:,1),AngleForceTorque(:,2));
xlabel('Turntable Angle (rad)')
ylabel('Force in x-direction (N)')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on


figure(3)
l=plot(AngleForceTorque(:,1),AngleForceTorque(:,7));
xlabel('Turntable Angle (rad)')
ylabel('Torque about z-direction (N m)')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on
