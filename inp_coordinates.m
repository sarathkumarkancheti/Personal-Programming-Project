%% Input coordinates
%__________________________________________________________________________
% INPUT: Number of coordinates, problem number, .INP files
%__________________________________________________________________________
% OUTPUT: Returns model parameters node numbers, element connectivity
% matrix, nodal coordinates, boundary nodes
%__________________________________________________________________________
function [nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
    right_edge_nodes,top_edge_nodes,nodeCoords,flux_nodes] = inp_coordinates(ncoord,problem)
% File path
if ncoord == 2 && problem == 1
    % for mesh with 4 node elements
    inpFile = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\1-D.inp";
elseif ncoord ==2 && problem == 4 || problem == 2
    %for 2 dimensional mesh
   inpFile = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\quad_2.inp";
elseif ncoord ==3 && problem == 3
    %for 3 dimensional mesh
   inpFile = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\3D_mesh_fine.inp"';
elseif ncoord == 2 && problem == 5
    % for mesh with two solvents
   inpFile = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\two_solvent_mesh.inp";
elseif ncoord == 2 && problem == 6
 % For patch test
 % Nodal coordinates
 nodeCoords = zeros(9,2);
 nodeCoords(1,:) = [0,0];
 nodeCoords(2,:) = [1,0];
 nodeCoords(3,:) = [2.5,0];
 nodeCoords(4,:) = [0,1];
 nodeCoords(5,:) = [1.2,1.2];
 nodeCoords(6,:) = [2.5,1];
 nodeCoords(7,:) = [0,2];
 nodeCoords(8,:) = [0.7,2];
 nodeCoords(9,:) = [2.5,2];
 % Element_connectivity
 elementNodes = zeros(4,4);
 elementNodes(1,:) = [1,2,5,4];
 elementNodes(2,:)= [2,3,6,5];
 elementNodes(3,:) = [4,5,8,7];
 elementNodes(4,:) = [5,6,9,8];
 nodeIDs = [1,2,3,4,5,6,7,8,9];
 elementIDs = [1,2,3,4];
 nnodes = 9;
 left_edge_nodes = [1,4,7]; bottom_edge_nodes = [1,2,3]; right_edge_nodes = [3,6,9]; top_edge_nodes = [7,8,9]; flux_nodes = [];
 flux_nodes = top_edge_nodes;
return;
    
end
% Open the file
fid = fopen(inpFile, 'r');
if fid == -1
    error('Could not open file');
end

% Initialize storage
nodes = [];    % Node coordinates
elements = []; % Element connectivity

% Flags for parsing
parsingNodes = false;
parsingElements = false;

% Read the file line by line
while ~feof(fid)
    line = strtrim(fgetl(fid)); 
    % Check for keywords
    if startsWith(line, '*Node', 'IgnoreCase', true)
        parsingNodes = true;
        parsingElements = false;
        continue;
    elseif startsWith(line, '*Element', 'IgnoreCase', true)
        parsingNodes = false;
        parsingElements = true;
        continue;
    elseif startsWith(line, '*Nset', 'IgnoreCase', true)
        parsingNodes = false;
        parsingElements = false;
        continue;
    elseif startsWith(line, '*', 'IgnoreCase', true)
        % End current section if another keyword starts
        parsingNodes = false;
        parsingElements = false;
        continue;
    end
    
    % Parse data
    if parsingNodes
        data = str2double(strsplit(line, ','));
        nodes = [nodes; data]; %#ok<AGROW>
    elseif parsingElements
        data = str2double(strsplit(line, ','));
        elements = [elements; data]; %#ok<AGROW>
    end
end

% Close the file
fclose(fid);

% Separate node data
nodeIDs = nodes(:, 1);
nodeCoords = nodes(:, 2:end);
nnodes = length(nodes(:,1));
% Separate element data
elementIDs = elements(:, 1);
elementNodes = elements(:, 2:end);
% Flux nodes for hydrogel block pumped with fluid 2D case
flux_nodes = nodes(nodeCoords(:, 1) <= 0.61 & nodeCoords(:, 2) == 1, 1);

if problem == 5
left_edge_nodes = nodes(nodeCoords(:,1)==0);
bottom_edge_nodes = nodes(nodeCoords(:,2)==0);
right_edge_nodes = nodes(nodeCoords(:,1)==1);
top_edge_nodes = nodes(nodeCoords(:,2)==2);
else
left_edge_nodes = nodes(nodeCoords(:,1)==0);
bottom_edge_nodes = nodes(nodeCoords(:,2)==0);
right_edge_nodes = nodes(nodeCoords(:,1)==1);
top_edge_nodes = nodes(nodeCoords(:,2)==1);
end
    
%{
disp('Node Coordinates:');
disp(nodeCoords);
disp('Element Connectivity:');
disp(elementNodes);
%}
end