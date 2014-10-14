function VB5_license(funcName)
% res=VB5_license(funcName)
%?
% Prints a short license text to the command line.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_license.m, prints a short license text for the vbSPT package
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén, bmelinden@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code

fprintf(...
['\nvbSPT, %s, Copyright (C) 2014 Martin Lindén\n' ...
 'v1 Martin Lindén and Fredrik Persson 2012.\n' ...
 'v2 with reparameterized transition matrix, Martin Lindén 2013.\n' ...
 'v3 with model of blur/localization errors, Martin Lindén 2014.\n\n' ...
 'This program comes with ABSOLUTELY NO WARRANTY. \n' ...
 'This is free software, and you are welcome to redistribute it \n' ...
 'under certain conditions. See license.txt for details. \n\n' ...
 'Additional permission under GNU GPL version 3 section 7 \n\n'...
 'If you modify this Program, or any covered work, by linking or combining it \n'...
 'with Matlab or any Matlab toolbox, the licensors of this Program grant you \n'...
 'additional permission to convey the resulting work. \n\n'],funcName);

end



