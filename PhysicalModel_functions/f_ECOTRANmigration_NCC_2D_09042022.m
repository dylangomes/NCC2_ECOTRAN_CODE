function ECOTRANmigration = f_ECOTRANmigration_NCC_2D_09042022(ECOTRANphysics)
% returns area of overlap (m2) between sub-domain boxes to which are later
% applied swimming speeds between boxes to get migration volumetric flux
% rates (m3/d).
% Overlaps are defined manually for each model.
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       ECOTRANphysics
%           num_t                   number of time steps	(vertical vector: num_t X 1)
%           num_boxes               1-MLD, 2-epipelagic, 3-mesopelagic, 4-bathypelagic, 5-benthic (QQQ delete MLD)
%           BoxArea_y               side or face area; North-South (y) orientation; (m2)	(2D matrix: time X num_boxes)
%           BoxArea_x               side or face area; East-West (x) orientation; (m2)      (2D matrix: time X num_boxes)
%           BoxFloorArea            floor area; (m2)        (2D matrix: time X num_boxes SOURCE)
%
% returns:
%       ECOTRANmigration
%           num_boxes
%           num_fluxes_migration	migration fluxes (potential number per group)
%           MIGRATION               migration as boundary area between boxes; (m2); (3D matrix: num_t X source (num_boxes+1) X destiny (num_boxes+1))
%           fname_ECOTRANmigration	name of this f_ECOTRANmigration function
%
% revision date: 11-25-2020
%       9/4/2022    changed name to specify that this function is for the 2D 5-box model case


% *************************************************************************
% STEP 1: basic domain parameters (GoMexOcn)-------------------------------
% step 1a: function code name ---------------------------------------------
fname_ECOTRANmigration	= mfilename; % name of this f_ECOTRANmigration function
display(['Running: ' fname_ECOTRANmigration])

% step 1b: unpack input terms ---------------------------------------------
num_t                   = ECOTRANphysics.num_t;         % number of time steps; (vertical vector: num_t X 1)
num_boxes               = ECOTRANphysics.num_boxes;     % (scalar: num_boxes); I, II, III, IV, V
% BoxArea_y               = ECOTRANphysics.BoxArea_y;	  % side or face area; North-South (y) orientation; (m2); (2D matrix: time X num_boxes)
% BoxArea_x               = ECOTRANphysics.BoxArea_x;	  % side or face area; East-West (x) orientation; (m2); (2D matrix: time X num_boxes)
BoxFloorArea            = ECOTRANphysics.BoxFloorArea;	% floor area; (m2); (2D matrix: num_t X SOURCE (num_boxes))
num_fluxes_migration	= 4;                            % SSS migration fluxes (potential number per group); QQQ only considering vertical migration at present
% *************************************************************************



% *************************************************************************
% STEP 2: Physical Structure (GoMexOCN)------------------------------------
%         define spatial relationships between boxes

% step 2a: vertical migration ---------------------------------------------
vertical_links          = [0 0 0 0 0 0
                           0 0 1 0 0 0
                           0 1 0 0 0 0
                           0 0 0 0 1 0
                           0 0 0 1 0 0
                           0 0 0 0 0 0]; % (2D matrix: SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))

vertical_links          = reshape(vertical_links', [1, (num_boxes+1), (num_boxes+1)]); % (3D matrix: 1 X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary)); NOTE transpose
vertical_links          = repmat(vertical_links, [num_t, 1, 1]); % (3D matrix: num_t X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))
vertical_links(:, 2, 3)	= BoxFloorArea(:, 2);
vertical_links(:, 3, 2)	= BoxFloorArea(:, 2);
vertical_links(:, 4, 5)	= BoxFloorArea(:, 4);
vertical_links(:, 5, 4)	= BoxFloorArea(:, 4);
% -------------------------------------------------------------------------


% step 2b: east-west migration --------------------------------------------

% step 2c: north-south migration ------------------------------------------

% step 2d: combine all migration links ------------------------------------
migration_links         = vertical_links; % add other links here (east-west face overlap & north-south face overlap); (3D matrix: num_t X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))
% *************************************************************************



% *************************************************************************
% STEP 3: define flux time-series------------------------------------------
MIGRATION               = migration_links;	% migration as boundary area between boxes; (m2); (3D matrix: num_t X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))
% *************************************************************************



% *************************************************************************
% STEP 4: pack results for export-----------------------------------------
ECOTRANmigration.num_boxes              	= num_boxes;
ECOTRANmigration.num_fluxes_migration       = num_fluxes_migration;     % migration fluxes (potential number per group)
ECOTRANmigration.MIGRATION                  = MIGRATION;                % migration as boundary area between boxes; (m2); (3D matrix: num_t X SOURCE (num_boxes+boundary) X DESTINY (num_boxes+boundary))
ECOTRANmigration.fname_ECOTRANmigration     = fname_ECOTRANmigration;	% name of this f_ECOTRANmigration function
% *************************************************************************


% end m-file***************************************************************