function freewakeSetup = defaultRun()
    omega = 176;
    freewakeSetup.integrationScheme = 'PCC';
    freewakeSetup.omega = 0.0; %Roughly mach 0.6 for a 7.5 ft span rotor (Caradonna Tung Rotor)
    freewakeSetup.dt = 2*pi/omega/50;
    freewakeSetup.nt = 650;
    freewakeSetup.globalLinearVelocity = [0, 0, 0];
    freewakeSetup.globalRotationAxis = [0, 0, 1];
    
    freewakeSetup.refL = 6.0;
    freewakeSetup.refC = 1.0;
    freewakeSetup.refA = pi*freewakeSetup.refL^2;
    freewakeSetup.refV = omega * freewakeSetup.refL;
    freewakeSetup.vMach = 1100;
    
    freewakeSetup.nSurfaces = 2;
    freewakeSetup.nChord = 9;
    freewakeSetup.nSpan = 22;
    freewakeSetup.nNearWake = 250; %170
    freewakeSetup.nFarWake = 2;
    
    freewakeSetup.surfaceAR = 5.25;
    freewakeSetup.surfacePitch = 5.0*pi/180;
    freewakeSetup.surfaceRotationRate = omega;
    freewakeSetup.surfaceTipDihedral = 15*pi/180;
    freewakeSetup.surfaceTipDihedralBreak = 0.90*5.25/6.0;
    
    freewakeSetup.isFreeWake = true;
    freewakeSetup.isFreeTipVortex = true;
    freewakeSetup.hasFixedTrailers = false;
    freewakeSetup.doPrandtlGlauert = true;
end