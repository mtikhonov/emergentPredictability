classdef (Abstract) V2communityFunctionLandscape < handle
    properties (SetAccess='protected', GetAccess='public')
        DIM     % dimension (number of species/mutations)
        cube  % storing the computed values of F. Sparse array with space for 2^12 entries
    end
    methods
        function obj = V2communityFunctionLandscape(DIM)
            if nargin<1
                DIM=10;
            end
            obj.DIM = DIM;
            if DIM>12 
                maxSize = 2^12;
            else
                maxSize = 2^DIM;
            end
            obj.cube = sparse([],[],[],2^DIM,1,maxSize);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = computeAll(obj,STOP_ON_FAILURE)
            assert(obj.DIM<=12,'computeAll is meant to be used with DIM <= 12 only.')
            if nargin<2
                STOP_ON_FAILURE = false;
            end
            C = dec2bin(0:(2^obj.DIM-1))=='1'; % all binary combinations
            f = getF(obj, C, STOP_ON_FAILURE);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = getF(obj,C,STOP_ON_FAILURE) % C = Community signature, binary vector
            if nargin<3
                STOP_ON_FAILURE = false;
            end
            assert(size(C,2)==obj.DIM,'Dimension mismatch')
            f = queryCube(obj, C);
            % find values not yet computed
            needToCompute = isnan(f);
            Csubset = C(needToCompute,:);
            NN = size(Csubset,1);
            valuesToCompute = f(needToCompute); % all NaNs
            parfor i=1:NN
            %for i=1:NN
                valuesToCompute(i) = getF_internal(obj,Csubset(i,:));
                if STOP_ON_FAILURE && isnan(valuesToCompute(i))
                    error('getF:fail','Stopping on failure as requested');
                end
            end
            storeInCube(obj, Csubset, valuesToCompute);
            f(needToCompute) = valuesToCompute;
        end
        
        % Helper functions:
        function f = queryCube(obj, C1)
            % C1 -> linear index
            ind = C2ind(obj, C1);
            % for a sparse array, absent values are 0
            % for us, 0 could actually be the value stored
            % whereas absent value means NaN
            f = full(obj.cube(ind));
            f = obj.swapNAN(f);
        end
        function ind = C2ind(obj, C)
            pow2 = 2.^uint32(0:obj.DIM-1);
            ind = (uint32(C))*pow2(:)+1;
        end
        function storeInCube(obj, C, f)
            f = obj.swapNAN(f);
            ind = C2ind(obj, C);
            obj.cube(ind) = f;
        end
        function [fPred] = leaveOneOut(obj)
            allC = dec2bin(0:(2^obj.DIM-1))=='1'; % all binary combinations
            fPred = NaN(size(1, 2^obj.DIM));
            for ii=1:length(fPred)
                % drop that point & train regression
                C1 = allC(ii,:);
                Cbut1 = allC;
                Cbut1(ii,:) = [];
                coef = trainRegression(obj, Cbut1);
                [~, ~, ~, fPred(ii)] = testRegression(obj, coef, C1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The regression training/testing functions take either the number
        % & density of the desired set as inputs (three-argument format):
        %     obj.trainLasso(N, p) 
        % or the actual set to use (two-argument format):
        %     obj.trainLasso(communitySet) 
        function [coef, CTrain, fTrain] = trainLasso(obj, N, p) 
            if nargin>=3
                CTrain = obj.randomSet(N,p);
            else
                CTrain = N;
            end
            fTrain = getF(obj,CTrain);
            [B, FitInfo] = lasso(double(CTrain), fTrain, 'Standardize', false, 'CV', 3);
            idxLambda1SE = FitInfo.Index1SE;
            coef = [FitInfo.Intercept(idxLambda1SE); B(:,idxLambda1SE)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [coef, CTrain, fTrain] = trainRegression(obj, N, p) % 2nd output = communities used for training
            if nargin>=3
                CTrain = obj.randomSet(N,p);
            else
                CTrain = N;
            end
            N = size(CTrain,1);
            fTrain = getF(obj,CTrain);
            coef = [ones(N,1), double(CTrain)]\fTrain;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [c, RMSE, CTest, fPredicted] = testRegression(obj, coef, N, p)
            if nargin>3
                CTest = obj.randomSet(N,p);
            else
                CTest = N;
            end
            N = size(CTest,1);
            fTest = getF(obj,CTest);
            fPredicted = [ones(N,1), double(CTest)] * coef(:);
            RMSE = sqrt(mean((fTest-fPredicted).^2));
            c = corr(fTest(:),fPredicted(:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [c, regressionCoef, subsetsArray] = noDiagsRegressionGridHeatmap(obj, N, pRange)
        % This function does not run regression on the diagonal of the
        % matrix. 
        %
        % Assay prediction power of regression trained at one compelxity
        % for making predictions at another. For this, assay a range of
        % training & testing sets of densities in pRange (a vector).
        
        % To reduce the # of communities that need to be equilibrated, select ONE
        % training set of each complexity. Use this same set for training a
        % regression model AND for evalauting other regression models. FOr the
        % diagonal entries, this would mean that the performance we compute is
        % for in-sample prediction. Which is okay: for out-of-sample prediction, we
        % can look at the pixels right next to the diagonal. And comparing the
        % diagonal pixels with neighboring ones will tell us how generalizable
        % the regression predictions are.
            % start by assembling the training/testing sets
            % and training the corresponding regression models
            NP = length(pRange);
            subsetsArray = cell(1,NP);
            regressionCoef = cell(1,NP);
            fprintf('Training regressions')
            for pp=1:NP
                fprintf('.');
                subs = obj.randomSet(N,pRange(pp));
                subsetsArray{pp} = subs;
                regressionCoef{pp} = obj.trainRegression(subs);
            end
            fprintf('\n');
            % Now, fill the grid testing in j the regression trained on i:
            c = NaN(NP, NP);
            for i=1:NP
                fprintf('.');
                for j=1:NP
                    c(i,j) = obj.testRegression(regressionCoef{i},subsetsArray{j});
                end
            end
        end

        function [c, regressionCoef, diag_testing_subsetsArray, ...
                training_subsetsArray] = regressionGridHeatmap(obj, N, pRange)
        % Same as the above regressionGridHeatmap but this time the
        % regression is tested on a unique set of communities on the
        % diagonal. As a result, it is slower. 
        %
        % N is the number of communities generated for training/testing
        % pRange is the size of communities being tested / total species
        %    This is the probability of a species being in a community that
        %    results in the desired community size. (If you want on average
        %    2 species in a community, then odds are 2/45 a species is in
        %    community. 
        %
        % Assay prediction power of regression trained at one compelxity
        % for making predictions at another. For this, assay a range of
        % training & testing sets of densities in pRange (a vector).
        
        % To reduce the # of communities that need to be equilibrated, select ONE
        % training set of each complexity. Use this same set for training a
        % regression model AND for evalauting other regression models, 
        % EXCEPT on the diagonals where a separate training set is created.

        % The expected number of communities tested in this one is double
        % the non-diagonal method. 
            
            % start by assembling the training/testing sets
            % and training the corresponding regression models
            NP = length(pRange); %Number of pixels on matrix
            %Contains info for each community size + special cell for
            %diagonals
            training_subsetsArray = cell(1,NP); 
            diag_testing_subsetsArray = cell(1,NP);
            %Will contain regression coefficients for each community size
            regressionCoef = cell(1,NP); 
            fprintf('Training regressions')
            for pp=1:NP %For each community size: 
                fprintf('.'); 
                % Generate testing/training communitites
                subs = obj.randomSet(N,pRange(pp));
                training_subsetsArray{pp} = subs;
                % Generate separate testing subset for diagonals
                subs = obj.randomSet(N,pRange(pp));
                % Remove communities used to train the diagonal
                [~,flagged] = intersect(subs, ...
                    training_subsetsArray{pp},'rows'); %Find matching rows
                subs(flagged,:) = []; %Removing matching rows
                diag_testing_subsetsArray{pp} = subs;
                
                %Train on the training communities
                regressionCoef{pp} = obj.trainRegression...
                    (training_subsetsArray{pp});
            end
            fprintf('\n');
            % Now, fill the grid testing in j the regression trained on i:
            c = NaN(NP, NP);
            disp('Now beginning testing.')
            for i=1:NP
                fprintf('.');
                for j=1:NP
                    % If on the diagonal, use diagonal testing set.
                    if j == i
                        c(i,j) = obj.testRegression...
                            (regressionCoef{i},diag_testing_subsetsArray{j});
                    % Else use one of the training sets.
                    else
                        c(i,j) = obj.testRegression...
                            (regressionCoef{i},training_subsetsArray{j});
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function C = randomSet(obj,N,p)
%             % This version picks N random subsets where each has exactly
%             % K=p*obj.DIM species - but such examples are rank-deficient
%             vals = rand(N, obj.DIM);
%             K = round(p*obj.DIM);
%             [~,speciesIdx] = sort(vals,2);
%             speciesIdx = speciesIdx(:,1:K);
%             subsetIdx = repmat(1:N, [K,1])';
%             C = false(N,obj.DIM);
%             C(sub2ind(size(C), subsetIdx(:), speciesIdx(:))) = true;
%         end
        function C = randomSet(obj,N,p)
            C = rand(N,obj.DIM)<p; %Random matrix with communities of approximately right size
            before = size(C,1);
            % remove duplicates
            C = unique(C,'rows'); %Removes duplicate communities
            % remove empty set 
            notNull = any(C,2); %Finds rows that aren't all zeros
            C = C(notNull, :); %Removes empty rows
            after = size(C,1);
            if after<before
                warning('Removed %d duplicate entries from regression training set', before-after);
            end
        end
    end
    methods (Static)
        function f = swapNAN(f)
            zers = f==0;
            nans = isnan(f);
            f(zers) = NaN;
            f(nans) = 0;
        end
    end
    
    methods (Abstract, Access='protected')
        getF_internal(obj,C1) % C1 = a single community signature, returns a number        
    end
end