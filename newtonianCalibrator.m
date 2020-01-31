%All units are SI, always.
% -------------- Configure Code Behaviour --------------------
FancyGraphics = 'False';
DisplayFigures = 'False';
SaveFigures = 'False';
SaveData = 'False';

% ----------- Load LHO DARM Spectrum -------------------------
%rawDARM=load('LHO_strain_noise.mat');
%FDARM=rawDARM.output(:,1);
%DARM=rawDARM.output(:,2);
Vec1F = [];
Vec2F = [];
Vec3F = [];
for j=1:20
    
    rotF=6; % Rotation frequency
    tTotal=100; % Integration time

    a2l=0.1; % Angle to Length coupling

    % NCal Cylinder Parameters
    CylinderHeight = 2*0.0254;
    CylinderDiameter = 1.5*0.0254;
    CylinderMass = 1.0558;
    CylinderAxialGridPoints = j;
    CylinderRadialGridPoints = j;

    % NCal Rotor Parameters
    RotorRadius2 = 2.375*0.0254; % Quadrople radius
    RotorRadius3 = 4.125*0.0254; % Hexapole radius

    RotorPosition = [-1.12365 -1.12365 0]-[-0.2 -0.39946 0];

    % Test Mass Parameters
    TMLength = 30e-2;
    TMDiameter = 40e-2;
    TMMass = 40;
    TMMoment = 1/12*TMMass*(3*TMDiameter^2/4+TMLength^2/4);
    TMAxialGridPoints = j;
    TMRadialGridPoints = j;


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

    for angle = 0:2*pi/100:2*pi

        if getappdata(f,'canceling')
            break
    end
    
	RotatedRotor = rotatePMArray(Rotor, angle, [0 0 1]);
	TranslatedRotatedRotor = translatePMArray(RotatedRotor, RotorPosition);    
    
    if FancyGraphics == 'True'
        displayPoints(TM, TranslatedRotatedRotor);
    end 
    
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

    tim=angle/2/pi/rotF;
    sampF=1/(tim(2)-tim(1));

   force=AngleForceTorque(:,2);%Force in x-direction
    torque=AngleForceTorque(:,7);%Torque about z-direction

    force=force-mean(force);
    torque=torque-mean(torque);
    
    TimeForce = [tim  force];

   
    
%% Fitting
    x1=[cos(2*pi*rotF*tim) sin(2*pi*rotF*tim)];
    x2=[cos(2*pi*2*rotF*tim) sin(2*pi*2*rotF*tim)];
    x3=[cos(2*pi*3*rotF*tim) sin(2*pi*3*rotF*tim)];


    % Linear least sqaure fitting(acclearation)
    w1 = inv(x1'*x1)*x1'*force;
    w2 = inv(x2'*x2)*x2'*force;
    w3 = inv(x3'*x3)*x3'*force;


    disp(['1F Amplitude: ' num2str(sqrt(w1(1)^2+w1(2)^2)) ' Phase: ' num2str(atan2(w1(2),w1(1))*180/pi) ' deg'])
    disp(['2F Amplitude: ' num2str(sqrt(w2(1)^2+w2(2)^2)) ' Phase: ' num2str(atan2(w2(2),w2(1))*180/pi) ' deg'])
    disp(['3F Amplitude: ' num2str(sqrt(w3(1)^2+w3(2)^2)) ' Phase: ' num2str(atan2(w3(2),w3(1))*180/pi) ' deg'])

   Vec1F = [Vec1F sqrt((w1(1))^2+(w1(2))^2)];
   Vec2F = [Vec2F sqrt((w2(1))^2+(w2(2))^2)];
   Vec3F = [Vec3F sqrt((w3(1))^2+(w3(2))^2)];
end
%% Save Force Time Series To file
fid = fopen('Force_Grid_Point_20.txt','wt');
for ii = 1:size(TimeForce,1)
    fprintf(fid,'%g\t',TimeForce(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
%%

if DisplayFigures == 'True'
    %Plot of Force Amplitude vs. GridPoint number
    figure(5)
    Plot1F = semilogy(1:length(Vec1F),Vec1F,1:length(Vec1F),Vec2F,1:length(Vec1F),Vec3F);
    legend('1F','2F','3F')
    xlabel('Grid Point Number')
    ylabel('Force Amplitude')
    %ylim([1e-17 3e-16])


    %Plot of Force Error Percentage vs. Grid point number
    figure(6)
    PlotEndPoints = plot(1:length(Vec1F),100*(Vec1F/ Vec1F(end)-1),1:length(Vec1F),100*(Vec2F / Vec2F(end)-1),1:length(Vec1F),100*(Vec3F / Vec3F(end)-1))
    legend('1F','2F','3F')
    xlabel('Grid Point Number')
    ylabel('Force Error %')




    figure(3)
    l=plot(angle,AngleForceTorque(:,2));
    xlabel('Turntable Angle (rad)')
    ylabel('Force in x-direction (N)')
    set(l,'LineWidth',1.5);
    set(gca,'FontSize',16);
    set(l,'MarkerSize',16);
    grid on
end 
