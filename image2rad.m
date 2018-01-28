% function [I_q,q,I_q_err,goodPixels,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = image2rad(experiment,runnum,photon_avg,D,goodPixels,Wavelength,x0,y0,tempdq)
function [I_q,I_q_adu,q,I_q_err,goodPixels,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = image2rad(experiment,runnum,photon_avg,photon_tot_adu,D,goodPixels,Wavelength,x0,y0,tempdq) %bms, changed input/output

%Provides additional experimental parameters to radialAverage.
%[I,q] = image2rad(photon_avg, dq)
%   photon_avg     : 388x185x32 photon map
%   dq              : width of q-bands (optional)
%   I               : Intensity as a function of q
%   q               : q (inverse angstoms)
%   
%   Intensity is given in units of average nubmer of photons per solid
%   angle per shot.
%
%   All of our experiment-specific parameters are set here. If you need to
%   change something, you shouldn't need to dig deeper than this file.

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

%% Set the area for each pixel
% (Note: Edge pixels are elongated by 2.5x in one dimension.) 
% From Henrik's email to Vale dated: 2013.01.09 13:12 (EST).
A=ones(388,185,32)*110*110;
A(194,:,:) = 110*110*2.5;
A(195,:,:) = 110*110*2.5;

%% Turn image into I(q)
% [I_q,q,I_q_err,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = radialAverage(experiment,runnum,photon_avg,Wavelength,D,x0,y0,dq,dphi,coors,A,goodPixels);
[I_q,I_q_adu,q,I_q_err,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = radialAverage(experiment,runnum,photon_avg,photon_tot_adu,Wavelength,D,x0,y0,dq,dphi,coors,A,goodPixels); %bms, changed input/output

