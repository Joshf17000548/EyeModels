function cornea = cornea( eye )
% Returns the cornea sub-field of an eye model structure
%
% Syntax:
%  cornea = human.cornea( eye )
%
% Description:
%   
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   cornea                - Structure.
%

%% Front corneal surface
% The Drasdo & Fowler model describes a corneal front surface that has a
% radius of curvature of 7.8 mm, and an eccentricity of 0.5. This code
% obtains the semi-radii for this ellipse, although it requires assuming
% something strange about D&F's definition of "eccentricity". Thanks to
% Giovanni Montesano of the University of London for his help in sorting
% this out. In the Table of the Drasdo & Fowler paper, values for
% "horizontal" and "vertical" radii in mm are given for the cornea. The use
% of these values, however, result in a cornea with optical power that is
% far too low.

%     syms a b 
%     Equation for the eccentricity (with a strange squaring of the denom)
%     eqn1 = 0 == (a^2-b^2)/a^2;
%     Equation for the radius of curvature
%     eqn2 = 7.8 == (b^2)/a;
%     sol = solve([eqn1, eqn2]);
%     radii = [eval(sol.a(2)), eval(sol.b(2)), eval(sol.b(2))];


% Corneal radii
radii = [8   8  8]; % looks like it should be radii
% radii = [15.6000   15.6  15.6 ];

% Create the quadric
S = quadric.scale(quadric.unitSphere,radii);

% We set the center of the cornea front surface ellipsoid so that the axial
% apex (prior to rotation) is at position [0, 0, 0]
S = quadric.translate(S,[-radii(1) 0 0]);

% Find the moster anterior point of this quadric surface
X = quadric.mostAnteriorPoint( S );

% Store these values
cornea.front.S = quadric.matrixToVec(S);
cornea.front.side = 1;
cornea.front.boundingBox=[-5 X(1) -9 9 -9 9];
cornea.front.center=[-radii(1) 0 0];


%% Store the combined corneal surfaces
cornea.S = cornea.front.S;
cornea.boundingBox = cornea.front.boundingBox;
cornea.side = 1;
cornea.mustIntersect = 1;
cornea.index = [];
cornea.label = {'cornea.front'};
cornea.plot.color = {[0.5 0.5 0.75]};


end