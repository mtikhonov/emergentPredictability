classdef V2cflLotkaVolterraPathogen < V2communityFunctionLandscape
    properties (SetAccess='protected')
        S % # of species in pool (minus the pathogen)
        inters % Interaction matrix, S+1 x S+1  
        growth % Growth rates, S+1
        K % Carrying capacities, S+1
        abd0 %Starting abundances S+1

        Teq %Maximum time of while loop (iterations x dt = time)
        
        % PATHOGEN PROPERTIES SHOULD BE CONTAINED IN LAST ROW OF ALL
        % PARAMETERS
    end

    methods 
        % This constructor builds our community & calculates abundance of
        % the pathogen. 
        function obj = V2cflLotkaVolterraPathogen(params)
            if nargin<1 %Error handling if no parameters provided
                error('Can''t construct a class instance with no params supplied.')
            end
            % params contain the info needed for LV random model
            % S = number of species in the pool minus the pathogen
            % inters = interaction matrix (S+1 x S+1)
            % growth = growth vector, growth rate for each species
            % K = carrying capacity vector, carry. cap. for each species
            % abd0 = initial abundances vector, initial abd. for each
            %           species
            % Teq = maximum time allowed equilibrium 
            % Epsi = threshold for community becoming "stable"
            
            

            % call parent constructor
            % obj is kinda like self in python, and by having this function
            % be the same name as the class, this function gets called when
            % the class is declared.
            obj@V2communityFunctionLandscape(params.S);
            
            %% Setting up initial data
            %Pass params to object
            obj.S = params.S;
            obj.Teq = params.Teq;
            obj.inters = params.inters;
            obj.K = params.K(:); 
            obj.growth = params.growth(:);
            obj.abd0 = params.abd0(:);
        end
    end
    methods (Access = 'protected')
        %This function evolves our community w/ LV dynamics 
        function abds = evolve_lv_model(obj,C1)
            % C1 describes species other than the pathogen. 
            % The pathogen is always added.
            select = [C1(:); true]; % a column vector of species present

            % Restrict the full LV system to the species that are present
            n0 = obj.abd0(select);
            RR = obj.growth(select);
            AA = obj.inters(select, select);
            KK = obj.K(select);
            dndt = @(t,n) RR./KK .* n .* (KK - n - AA*n);
            [~, n] = ode15s(dndt, [0 obj.Teq], n0);
            
            abds = zeros(obj.S+1,1);
            abds(select) = n(end,:);
        end

        function f = getF_internal(obj,C1) 
            % C1 = a single binary community signature, returns a number
            abd = evolve_lv_model(obj, C1);
            f = abd(end); %Returns only the pathogen's final abundance
        end
    end
end