%All units are SI, always.

rawDARM=load('LHO_strain_noise.mat');
FDARM=rawDARM.output(:,1);
DARM=rawDARM.output(:,2);

rotF=30; % Rotation frequency
tTotal=100; % Integration time

a2l=0.1; % Angle to Length coupling

% NCal Cylinder Parameters
CylinderHeight = 2*0.0254;
CylinderDiameter = 1.5*0.0254;
CylinderMass = 1.0558;
CylinderAxialGridPoints = 5;
CylinderRadialGridPoints = 5;

% NCal Rotor Parameters
RotorRadius2 = 2.375*0.0254; % Quadrople radius
RotorRadius3 = 4.125*0.0254; % Hexapole radius

RotorPosition = [-1.12365 -1.12365 0]-[-0.2 -0.39946 0];

% Test Mass Parameters
TMLength = 30e-2;
TMDiameter = 40e-2;
TMMass = 40;
TMMoment = 1/12*TMMass*(3*TMDiameter^2/4+TMLength^2/4);
TMAxialGridPoints = 10;
TMRadialGridPoints = 10;


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

%Make a Test Mass
TM = genPointMassAnnlSheet(TMMass, 0, TMDiameter/2, TMLength,...
		TMAxialGridPoints, TMRadialGridPoints);

%Spin the rotor and compute the torques
AngleForceTorque = [];

f = waitbar(0,'Calculating','Name','NCal PointGravity',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for angle = 0:2*pi/1000:2*pi

    if getappdata(f,'canceling')
        break
    end
    
	RotatedRotor = rotatePMArray(Rotor, angle, [0 0 1]);
	TranslatedRotatedRotor = translatePMArray(RotatedRotor, RotorPosition);    

	displayPoints(TM, TranslatedRotatedRotor); 
	[Force Torque ] = pointMatrixGravity(TM, TranslatedRotatedRotor);

	AngleForceTorque = [AngleForceTorque; angle Force Torque];
    
    waitbar(angle/(2*pi),f);
end
delete(f);
%% Force and Torque vs Rotor Angle

Ncycles=tTotal*rotF;
angle=[];

for j=0:Ncycles-1
    
    angle=[angle; AngleForceTorque(:,1)+2*pi*j];
    
end

AngleForceTorque=repmat(AngleForceTorque,Ncycles);


%% Strain Calculations

tim=angle/2/pi/rotF;
sampF=1/(tim(2)-tim(1));

accel=AngleForceTorque(:,2)/TMMass;
angAccel=AngleForceTorque(:,7)/TMMoment;

accel=accel-mean(accel);
angAccel=angAccel-mean(angAccel);

s=0;
sAng=0;
vel=[];
velAng=[];

for j=1:length(accel)
    
    vel=[vel s+accel(j)*(tim(2)-tim(1))];
    s=s+accel(j)*(tim(2)-tim(1));
    
    velAng=[velAng sAng+angAccel(j)*(tim(2)-tim(1))];
    sAng=sAng+angAccel(j)*(tim(2)-tim(1));
end

vel=vel-mean(vel);
velAng=velAng-mean(velAng);

s=0;
sAng=0;
x=[];
xAng=[];

for j=1:length(vel)
    
    x=[x s+vel(j)*(tim(2)-tim(1))];
    s=s+vel(j)*(tim(2)-tim(1));
    
    xAng=[xAng sAng+velAng(j)*(tim(2)-tim(1))];
    sAng=sAng+velAng(j)*(tim(2)-tim(1));
end

x=x-mean(x);
xAng=xAng-mean(xAng);

strain=(x+a2l*xAng)/4e3;

%% Strain Spectrum Calculations
[AAng, ~] = asd2(accel,1/sampF, 1, 1, @hann);
[AX, F] = asd2(angAccel,1/sampF, 1, 1, @hann);

[AS, F2] = asd2(strain,1/sampF, 1, 1, @hann);

AStrain=(AX+a2l*AAng)./F'.^2/4/pi^2/4e3;

%% Fitting
x1=[cos(2*pi*rotF*tim) sin(2*pi*rotF*tim)];
x2=[cos(2*pi*2*rotF*tim) sin(2*pi*2*rotF*tim)];
x3=[cos(2*pi*3*rotF*tim) sin(2*pi*3*rotF*tim)];

% Linear least square fitting
w1=inv(x1'*x1)*x1'*strain';
w2=inv(x2'*x2)*x2'*strain';
w3=inv(x3'*x3)*x3'*strain';

disp(['1F Amplitude: ' num2str(sqrt(w1(1)^2+w1(2)^2)) ' Phase: ' num2str(atan2(w1(2),w1(1))*180/pi) ' deg'])
disp(['2F Amplitude: ' num2str(sqrt(w2(1)^2+w2(2)^2)) ' Phase: ' num2str(atan2(w2(2),w2(1))*180/pi) ' deg'])
disp(['3F Amplitude: ' num2str(sqrt(w3(1)^2+w3(2)^2)) ' Phase: ' num2str(atan2(w3(2),w3(1))*180/pi) ' deg'])
%%
figure(2)
l=loglog(F,AStrain,F2,AS,FDARM,DARM);
text(rotF,sqrt(w1(1)^2+w1(2)^2),' 1F')
text(2*rotF,sqrt(w2(1)^2+w2(2)^2),' 2F')
text(3*rotF,sqrt(w3(1)^2+w3(2)^2),' 3F')
xlabel('Frequency (Hz)')
ylabel('Strain (1/\surd(Hz))')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
xlim([1e1 1e3])
ylim([1e-24 1e-20])
grid on

figure(3)
l=plot(angle,AngleForceTorque(:,2));
xlabel('Turntable Angle (rad)')
ylabel('Force in x-direction (N)')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on

figure(4)
l=plot(angle,AngleForceTorque(:,7));
xlabel('Turntable Angle (rad)')
ylabel('Torque about z-direction (N m)')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on
