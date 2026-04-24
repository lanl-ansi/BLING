function [x, x0, par] = rescaling(x, x0, par)
 
% Pressures at nodes
x(par.index.nodes) = x(par.index.nodes) * par.varp;
x0(par.index.nodes) = x0(par.index.nodes) * par.varp;
par.xmin(par.index.nodes) = par.xmin(par.index.nodes) * par.varp;
par.xmax(par.index.nodes) = par.xmax(par.index.nodes) * par.varp;

% Flows on comps
x(par.index.comps) = x(par.index.comps) * par.varq;
x0(par.index.comps) = x0(par.index.comps) * par.varq;
par.xmin(par.index.comps) = par.xmin(par.index.comps) * par.varq;
par.xmax(par.index.comps) = par.xmax(par.index.comps) * par.varq;

% Flows on pipes
x(par.index.pipes) = x(par.index.pipes) * par.varq;
x0(par.index.pipes) = x0(par.index.pipes) * par.varq;
par.xmin(par.index.pipes) = par.xmin(par.index.pipes) * par.varq;
par.xmax(par.index.pipes) = par.xmax(par.index.pipes) * par.varq;

% Flows at gnodes
par.q = par.q * par.varq;
x(par.index.gnode_d) = x(par.index.gnode_d) * par.varq;
x(par.index.gnode_s) = x(par.index.gnode_s) * par.varq;
x0(par.index.gnode_d) = x0(par.index.gnode_d) * par.varq;
x0(par.index.gnode_s) = x0(par.index.gnode_s) * par.varq;

par.xmax(par.index.gnode_d) = par.xmax(par.index.gnode_d) * par.varq;
par.xmax(par.index.gnode_s) = par.xmax(par.index.gnode_s) * par.varq;

par.xmin(par.index.gnode_d) = par.xmin(par.index.gnode_d) * par.varq;
par.xmin(par.index.gnode_s) = par.xmin(par.index.gnode_s) * par.varq;
end

