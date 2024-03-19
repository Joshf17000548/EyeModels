function retina = retina(eye)
% Returns the retina sub-field of an eye model structure
%
% Syntax:
%  retina = human.retina( eye )
%
% Description:
%   The retinal surface (i.e., the vitreous chamber / posterior segment) is
%   modeled as an ellipsoidal quadric surface, with adjustments in size for
%   the spherical ametropia of the eye. Atchison 2006 provides the semi-
%   radii of a triaxial ellipsoid model of the retinal surface, with these
%   values varying by spherical ametropia:
%
%       Atchison, David A. "Optical models for human myopic eyes." Vision
%       research 46.14 (2006): 2236-2250.
%
%   The Atchison model is arranged around the visual axis, with the
%   vitreous chamber tilted and shifted. In the current model, the retinal
%   ellipsoid is aligned with and centered on the optical axis.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   retina                - Structure.
%
% Examples:



% An essentially spherical retina. I add a tiny bit of asphericity to allow
% the ellipsoidal coordinate system to work with this model.
retinaRadiiEmmetrope = [10.148 11.455 11.365];
retinaRadiiAmetropiaSlope = [-0.163 -0.043 -0.090];

% axialLength = 23.58 - 0.299*eye.meta.sphericalAmetropia;
sphericalAmetropia = eye.R(:,15)/-0.299 + 78.86;

retinaRadii = ...
    retinaRadiiEmmetrope + retinaRadiiAmetropiaSlope.* sphericalAmetropia;

retinaCenter = ...
    [-(eye.R(:,15) - retinaRadii(1)) 0 0];

retina.radii = retinaRadii;

% Assemble the components
S = quadric.scale(quadric.unitSphere,retinaRadii);
S = quadric.translate(S,retinaCenter);

% Assemble the system
retina.S = quadric.matrixToVec(S);
retina.side = -1;
retina.boundingBox = [-eye.R(:,15) -9 -13 13 -12 12];
retina.mustIntersect = 1;
retina.label = {'retina'};
retina.plot.color = {[0.75, 0.75, 0.75]};

end

