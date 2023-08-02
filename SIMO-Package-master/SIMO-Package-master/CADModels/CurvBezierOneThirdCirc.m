% -------------------------------------------------------------------------
% a one-third circle modelled by a parameterized Bezier curve
% -------------------------------------------------------------------------

%{
Copyright (C) <2014-2016>  <Khanh Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

clear
close all
clc

radius = 1; % radius of the circle
alfa = 1/3 * pi; % 1/2 segment angle

% control points
CtrlPts = zeros(4, 3);
CtrlPts(1 : 2, 1) = [-sin(alfa); cos(alfa)];
CtrlPts(1 : 2, 2) = [0; 1/cos(alfa)];
CtrlPts(1 : 2, 3) = [sin(alfa); cos(alfa)];

CtrlPts(1 : 2, :) = CtrlPts(1 : 2, :) * radius;
% weights
CtrlPts(4, :) = 1;
fac = cos(alfa);
CtrlPts(:, 2) = CtrlPts(:, 2) * fac;

p = size(CtrlPts, 2) - 1; % degree of the curve

% evalutate the bernstein basis functions
xi = linspace(0, 1, 101); % allocate parametric points
B = AllBernstein(p, xi);

% contruct the rational bezier curve in 4D space
Cw = CtrlPts * B';

% project this curve to 3D space
w = Cw(4, :);
C = bsxfun(@rdivide, Cw, w);

% plot basis functions
figure
set(gcf,'color','white')
axis equal
daspect([1 1 1])
set(gcf, 'color', 'white');
plot(xi, B, 'LineWidth', 1.5)

% plot the curve, control points and control polygon
figure
hold on
set(gcf, 'color', 'white')
% plot the curve
plot(C(1, :), C(2, :), 'b-', 'LineWidth', 1.5)
% plot control points
weights = CtrlPts(4, :);
plot(CtrlPts(1, :)./weights, CtrlPts(2, :)./weights, 'r.', 'MarkerSize', 20)
% plot control net
plot(CtrlPts(1, :)./weights, CtrlPts(2, :)./weights, 'k--')

axis off
axis equal





