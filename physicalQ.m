function [l,q,qerr,di,avg_angle,distances] = physicalQ(experiment,radii,D,L)
%Uses detector geometry to convert physical space to q space.
%[l,q,qerr] = physicalQ(radii,D)
%
%   radii: radial distances from center of detector [um]
%   D    : distance from front of Be window to face of detector [um]
%   l    : path length seen by detector [um]
%   q    : position in q-space
%   qerr : error in q; s = q +/- qerr
%   di   : effective distance to detector [um] : distance from center of
%          the scattering region to the detector's face
%
% Revised last: Early October, 2012, JMB
% Converted to a function: 2013.01.11, VCS

%Parameters
k0_keV = 1.23984187/L;
h_bar = 6.58211928e-16;
c = 2.99792458*10^18; %given in angstrom per second

%List of boundary objects, x,y in mm
%High limits
upstream_200u_aperture = [-179.48,.1];
downstream_200u_aperture = [-179.65,.1];
if strcmp(experiment,'L560');
    flange_drill_hole = [-189.84,3.175]; % 3.175mm for bore hole.
else
    flange_drill_hole = [-189.84,(3.175/2)]; % 3.175/2mm for bore hole if washer is jammed in.
end
flange_bore = [-200.01,19.18];
be_holder_upstream_limit = [-201.28,21.39];
be_holder_downstream_limit = [-203.28,25.39];

%%%%% EDITED HERE %%%%%
%cspad = [-249.33,93.56]; %This value needs to be verified.
cspad_x = be_holder_upstream_limit(1)-(D/1000); %gets the absolute position of the CSPAD using to Be window

%Low Limits
screw_tip_1 = [-193.19,1.04];
screw_tip_2 = [-193.57,1.42];
nut = [-198.2,2.76];
washer = [-199.66,3.16+.6];
%Master lists of these limits:
potential_upper_limits = [upstream_200u_aperture;downstream_200u_aperture;flange_drill_hole;flange_bore;be_holder_upstream_limit;be_holder_downstream_limit];
potential_lower_limits = [screw_tip_1;screw_tip_2;nut;washer];

%%%%% EDITED HERE %%%%%
%r_active_region = [10.77:.1:93.56]; %This is the portion of the CS-PAD that will detect scattering.
r_active_region = radii'/1000;

%This formats r_active_region values to be useful to find limiting objects.
r_values = r_active_region';

%Here the HIGH limiting points are determined:
%High limit comparison takes the difference between the points, making 6
%pairs of columns: a pair for each of the objects. Instead of solving for
%the tangent of the angle, it's sufficient to compare the  of the
%lengths to find the steepest. 

high_limit_comparison=zeros([numel(r_values) numel(potential_upper_limits)]);

for i=1:(numel(potential_upper_limits)/2) % X and Y
    high_limit_comparison(:,2*i-1) = ones(numel(r_values),1)*abs(cspad_x)-abs(potential_upper_limits(i,1)); % relative X
    high_limit_comparison(:,2*i) = r_values-potential_upper_limits(i,2); % relative Y
    high_limit_comparison(high_limit_comparison<0)=0;
end

all_far_ratios=zeros([numel(high_limit_comparison(:,1)) numel(potential_upper_limits)/2]);
for i=1:(numel(potential_upper_limits)/2)
    all_far_ratios(:,i)=high_limit_comparison(:,2*i)./high_limit_comparison(:,2*i-1); % Y/X
end

limiting_far_ratio=max(all_far_ratios,[],2);
far_angle=atand(limiting_far_ratio);

%Here the LOW limiting points are determined:
%Low limit comparison takes the difference between the points, making 6
%pairs of columns: a pair for each of the objects. Instead of solving for
%the tangent of the angle, it's sufficient to compare the ratio of the
%lengths to find the steepest. 

low_limit_comparison=zeros([numel(r_values) numel(potential_lower_limits)]);

for i=1:(numel(potential_lower_limits)/2) % X and Y
    low_limit_comparison(:,2*i-1) = ones(numel(r_values),1)*abs(cspad_x)-abs(potential_lower_limits(i,1)); % relative X
    low_limit_comparison(:,2*i) = r_values-potential_lower_limits(i,2); % relative Y
    for f=1:numel(low_limit_comparison(:,2*i))
        if low_limit_comparison(f,2*i) < 0
            low_limit_comparison(f,(2:2:8)) = 0;
        end
    end
end

all_close_ratios=zeros([numel(low_limit_comparison(:,1)) numel(potential_lower_limits)/2]);
for i=1:(numel(potential_lower_limits)/2)
    all_close_ratios(:,i)=low_limit_comparison(:,2*i)./low_limit_comparison(:,2*i-1); % Y/X
end

limiting_close_ratio=min(all_close_ratios,[],2);
close_angle=atand(limiting_close_ratio);

%Note: These angles are measured starting from the point on the detector
%and are therefore different from the scattering angle!
max_visible = r_values./tand(far_angle);
min_visible = r_values./tand(close_angle);

%upper_angle = 90-max_angle;
s_max = 2* ((k0_keV*1000)/(h_bar*c)) * sind(far_angle/2);

%lower_angle = 90-min_angle;
s_min = 2* ((k0_keV*1000)/(h_bar*c)) * sind(close_angle/2);

reaction_length = max_visible - min_visible;
s_range = s_min - s_max;
s_range(s_range<0)=0;
s_avg = (s_min+s_max)/2;


%%%%% EDITED HERE %%%%%
% construct returned variables
l   = reaction_length;
q   = s_avg;
qerr= s_range/2;
di  = ((max_visible + min_visible) / 2);
avg_angle = atand(r_values./di); % Refers to scattering angle.

% filter out nonsense values
l   (l    < 0) = 0;
q   (q    < 0) = 0;
qerr(qerr < 0) = 0;
di((l==0) | (q==0) | (qerr==0)) = 0;

% mm->um conversion
l  = l * 1000;
di = di * 1000;

% For focusing parameters:
focustocspad=screw_tip_1(1)-cspad_x;
distances.distfromfocus = -(di-1000*(focustocspad));
min(distances.distfromfocus)
max(distances.distfromfocus)
distances.distfromfocus(distances.distfromfocus>0) = 0;

distances.distintochamber = 1000*(upstream_200u_aperture(1)-cspad_x)-di;
distances.distintochamber(distances.distintochamber<0) = 0;

%This part plots the reaction length:
figure(101);
scatter(s_avg(:),l(:))
xlabel('S, A^-^1')
ylabel('Length, mm')
legend('Reaction length')
%axis([0 100 0 10])

%The above plots the full graph with max, min, and range of s.
figure(102);
plot(r_active_region,s_range,r_active_region,s_max,r_active_region,s_min,r_active_region,q)
xlabel('Distance from Detector Center, mm')
ylabel('S-Values visible, A^-^1')
legend('s range','lower s limit','upper s limit','average s')
%plot(r_lower_region,upper_angle_result)

figure(103);
plot(r_active_region,2*qerr)
xlabel('Distance from Detector Center, mm')
ylabel('S-Value range visible, A^-^1')
legend('s range')
%axis([0 100 0 .7])
end
