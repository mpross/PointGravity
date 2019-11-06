% #    PointGravity: a simple point-wise Newtonian gravitation calculator.
% #    Copyright (C) 2017  Charles A. Hagedorn
% #
% #    This program is free software: you can redistribute it and/or modify
% #    it under the terms of the GNU General Public License as published by
% #    the Free Software Foundation, version 3.
% #
% #    This program is distributed in the hope that it will be useful,
% #    but WITHOUT ANY WARRANTY; without even the implied warranty of
% #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% #    GNU General Public License for more details.
% #
% #    You should have received a copy of the GNU General Public License
% #    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%Companion function to pointMassGravity
%Displays test mass (array1) in blue, field mass (array2) in red
function displayPoints(array1, array2)

    figure(1)
	plot3(0,0,0,'x', ...
        array1(:,2) , array1(:,3) , array1(:,4) , '.' ,...
	      	array2(:,2) , array2(:,3) , array2(:,4) , '.' );
	xlabel( 'x (m)' );
	ylabel( 'y (m)' );
	zlabel( 'z (m)' );
    axis('equal')
    xlim([-1.5 0.5])
    zlim([-0.5 0.5])
    ylim([-1.5 0.5])
    grid on

end
