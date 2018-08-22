function [R2,S,setup, Sinc, SIG, ir, R16] = calcR2(H,T,slope,igflag);
%
% [R2,S,setup, Sinc, SIG, ir] = calcR2(H,T,slope,igflag);
%
% Calculated 2% runup (R2), swash (S), setup (setup), incident swash (Sinc)
% and infragravity swash (SIG) elevations based on parameterizations from runup paper
% also Iribarren (ir)
% August 2010 - Included 15% runup (R16) statistic that, for a Guassian distribution, 
% represents mean+sigma. It is calculated as R16 = setup + swash/4.  
% In a wave tank, Palmsten et al (2010) found this statistic represented initiation of dune erosion. 
%
%
% H = significant wave height, reverse shoaled to deep water
% T = deep-water peak wave period
% slope = radians
% igflag = 0 (default)use full equation for all data
%        = 1  use dissipative-specific calculations when dissipative conditions exist (Iribarren < 0.3)
%        = 2  use dissipative-specific (IG energy) calculation for all data
%
% based on:
%  Stockdon, H. F., R. A. Holman, P. A. Howd, and J. Sallenger A. H. (2006),
%    Empirical parameterization of setup, swash, and runup,
%    Coastal Engineering, 53, 573-588.
% author: hstockdon@usgs.gov

% fix up the inputs
if ~exist('igflag','var')
    igflag = 0; 
end
H = H(:)      
T = T(:); 
slope =  slope (:) 
N = length(H);
g = 9.81;

% intialize output
R2 = nan*ones(N,1);
R16 = nan*ones(N,1);
S = nan*ones(N,1);
setup = nan*ones(N,1);
Sinc = nan*ones(N,1);
SIG = nan*ones(N,1);

% make slopes positive!
slope = abs(slope);

% compute wavelength and Iribarren
L = (g*T.^2) / (2*pi);
ir = slope./sqrt(H./L);

% do it
for ii = 1:N  % have to do in loop to check Iribarren number
    if igflag == 2      % use dissipative equations (IG) for ALL data
        R2(ii) = 1.1*( 0.039 * sqrt(H(ii).*L(ii)));
        S(ii) = 0.046*sqrt(H(ii).*L(ii));
        setup(ii) = 0.016*sqrt(H(ii).*L(ii));
        
    elseif igflag== 1 & ir(ii)<0.3		% if dissipative site use diss equations
        R2(ii) = 1.1*( 0.039 * sqrt(H(ii).*L(ii)));
        S(ii) = 0.046*sqrt(H(ii).*L(ii));
        setup(ii) = 0.016*sqrt(H(ii).*L(ii));
        
    else 	% if int/ref site, use full equations
        %  Coded as written in paper (equation 19)
        %         part1 = 0.35*slope(ii).*sqrt(H(ii).*L(ii));
        %         part2 = sqrt(H(ii).*L(ii).*( (0.5625*slope(ii).^2)+0.0036))/2;
        %         R2(ii) = 1.1*(part1+part2);
        setup(ii) = 0.35*slope(ii).*sqrt(H(ii).*L(ii));
        Sinc(ii) = 0.75*slope(ii).*sqrt(H(ii).*L(ii));
        SIG(ii) = 0.06*sqrt(H(ii).*L(ii));
        S(ii) = sqrt(Sinc(ii).^2 + SIG(ii).^2);
        R2(ii) = 1.1*(setup(ii) + S(ii)/2)
        R16(ii) = 1.1*(setup(ii) + S(ii)/4)
    end
end
