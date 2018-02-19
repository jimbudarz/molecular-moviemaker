function [I_q,I_q_adu,q,I_q_err,goodPixels,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = image2rad(experiment,runnum,photon_avg,photon_tot_adu,D,goodPixels,Wavelength,x0,y0,tempdq)
% Wrapper function for radialAverage. This provides experiment-specific
% parameters to describe the geometry of the experiment. All of our
% experiment-specific parameters are set here. If you need to change 
% something, you shouldn't need to dig deeper than this file.
%
%       experiment : The proposal number for the experiment (e.g. 'L560','560','l560','56012')
%           runnum : The run (a.k.a. scan) number to be loaded and processed 
%       photon_avg : 388x185x32 image matrix of photons on CSPAD
%   photon_tot_adu : 388x185x32 image matrix of total counts in ADU on CSPAD
%                D : Distance from the beryllium window to the CSPAD surface
%       goodPixels : 388x185x32 boolean matrix of pixels that are reliable
%       Wavelength : The wavelength of the incident X-ray in nanometers
%               x0 : Horizontal position of the x-ray beam on the detector (center of radial average)
%               y0 : Vertical position of the x-ray beam on the detector (center of radial average)
%           tempdq : width of q-bins (optional)
%              I_q : Intensity in photons per solid angle per shot as a function of q
%          I_q_adu : Intensity in ADU (analog-digital units) as a function of q
%                q : q (inverse angstoms)
%          I_q_err : Error in the intensity
%            qSort : Cell array of the pixel numbers that lie within each bin of q
%            I_phi : Intensity as a function of phi (azimuthal angle)
%        phivalues : Array of phi values in degrees
%    pctExcitation : The percent of excited molecules as a function of q.
%         angleMap : 388x185x32 matrix of average angle (theta) exposed to each pixel
%           phiMap : 388x185x32 matrix of average phi exposed to each pixel
%   

if strcmp(experiment,'b0114')
    coors=load('cspad_coordinates_b0114.mat');
    dq = 0.050; % inverse angstroms
    dphi = 2; % In degrees
%     x0 = 94370; % microns   (from findCSPADcenter.m, run 203)
%     y0 = 92520; % microns   (from findCSPADcenter.m, run 203)
elseif strcmp(experiment,'56012')
    coors=load('cspad_coordinates_L560.mat');
    dq = 0.050; % inverse angstroms
    dphi = 2; % In degrees
%     x0 = 94110; % microns   (Daniel's center x, L560)
%     y0 = 93960; % microns   (Daniel's center y, L560)
elseif strcmp(experiment,'i0613')
    coors=load('cspad_coordinates_i0613.mat');
    dq = 0.010; % inverse angstroms
    dphi = 2; % In degrees
%     x0 = 94110; % microns   (Daniel's center x, L560)
%     y0 = 93960; % microns   (Daniel's center y, L560)
end

if exist('tempdq','var')
    dq=tempdq;
end

%% Set the surface area for each pixel
% (Note: Edge pixels are elongated by 2.5x in one dimension.) 
% From Henrik's email to Vale dated: 2013.01.09 13:12 (EST).
A=ones(388,185,32)*110*110;
A(194,:,:) = 110*110*2.5;
A(195,:,:) = 110*110*2.5;

%% Turn image into I(q)
[I_q,I_q_adu,q,I_q_err,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = radialAverage(experiment,runnum,photon_avg,photon_tot_adu,Wavelength,D,x0,y0,dq,dphi,coors,A,goodPixels);

