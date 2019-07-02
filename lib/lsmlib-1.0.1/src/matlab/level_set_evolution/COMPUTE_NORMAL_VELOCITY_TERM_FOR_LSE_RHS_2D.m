%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COMPUTE_NORMAL_VELOCITY_TERM_FOR_LSE_RHS_2D() computes the 
% contribution of a normal velocity term to the right-hand side of 
% the level set equation.
%
% Usage:  lse_rhs = COMPUTE_NORMAL_VELOCITY_TERM_FOR_LSE_RHS_2D( ...
%                     normal_velocity, ...
%                     phi_x_plus, phi_y_plus, ...
%                     phi_x_minus, phi_y_minus)
%
% Arguments:
% - phi:               level set function
% - ghostcell_width:   ghostcell width for phi
% - normal_velocity:   normal velocity
% - phi_x_plus:        x-component of plus HJ ENO derivative
% - phi_y_plus:        y-component of plus HJ ENO derivative
% - phi_x_minus:       x-component of minus HJ ENO derivative
% - phi_y_minus:       y-component of minus HJ ENO derivative
%
% Return value:
% - lse_rhs:           normal velocity contribution to right-hand side 
%                      level set evolution equation
%
% NOTES:
% - The phi_x_plus, phi_y_plus, phi_x_minus and phi_y_minus arrays are 
%   assumed to be the same size
%
% - All data arrays are assumed to be in the order generated by the 
%   MATLAB meshgrid() function.  That is, data corresponding to the 
%   point (x_i,y_j) is stored at index (j,i).  The spatial derivative
%   functions provided by LSMLIB return data arrays ordered in this
%   way.
%
% - The returned lse_rhs array is the same size as phi.  However, only
%   the values of the RHS of the level set evolution equation within
%   the _interior_ of the computational grid are computed.  In other
%   words, values of the RHS in the ghostcells are _not_ computed; the
%   value in the ghostcells is set to 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyrights: (c) 2005 The Trustees of Princeton University and Board of
%                 Regents of the University of Texas.  All rights reserved.
%             (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:   $Revision: 149 $
% Modified:   $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

