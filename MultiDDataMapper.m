%> @file MultiDDataMapper.m
%>
%> @brief Simple class for creating a hashtable. Or at least that was 
%> the idea. An actual hashtable would just map any key to any value
%> while i am create multi dimensional matrices with multiple keys
%>
%> @section matlabComments Details
%
%> @authors Eoin O'Keeffe (eoin.okeeffe.09@ucl.ac.uk)
%> @date initiated: 17/08/2011
%> @version 
%> 1.0: Initially called HashTable. Now changed to MultiDDataMapper
%> <br />Version 2.0: Extended so it makes a bit more sense and there is
%> added functionality. However, it will remain backwards compatible
%>
%> @section intro Method
%> Takes in a multivariable array, and indexes it so it can be referenced
%> properly. Created for IEA for the large datasets. It has since been
%> extended to be used with GloTraM
%> @subsection version_history
%> 1.0 Initial set up for IEA
%> 2.0 Added to GloTraM and some commenting,intuitive code added. Remains
%> backwards compatible but data is now stored in table format
%>
%>@attention
%>@todo Finish getFilteredDataset so it aggregates over keys - can't find a
%> reason why i'd need it for now
classdef MultiDDataMapper < handle
    properties
        %> This is for backwards compatibility
        indxs
        keys
        
        %structure with one internal variable - keys
        keyStore
        %structure with one internal variable - indxs
        indxStore
        %> Stores the actual data in a multi dimensional matrix
        vals
        %> The number of keys
        noKeys
        %> The number of values
        noVals
        %> cell array of key names
        keyNames
        %> cell array of value names
        valNames
        %> stores the cardinality of the keys
        cardinality
        %> stores every in rows with each column corresponding to a
        %>dimension
        data
    end
    
    methods
        % ======================================================================
    %> @brief old function to set up the data
    %>
    %> @param obj
    %> @param data the actual data as a table
    %> @param 
    %> @retval
    % ======================================================================
        function setup5d(obj,data,values)
            warning('This will be removed in future versions. Please use the new  format');
           obj.indxs = zeros(size(data,1),size(data,2));
           obj.keys = obj.indxs;
           %Build referencing
           for i=1:size(data,2)
               tmp=sort(unique(data(:,i)),'ascend');
               obj.keys(1:size(tmp,1),i)=tmp;
               obj.indxs(1:size(tmp,1),i) = [1:size(tmp,1)]';
           end %for i
           obj.indxs = obj.indxs(any(obj.indxs,2),:); %remove rows that contain all 0s
           obj.keys = obj.keys(any(obj.keys,2),:); %remove rows that contain all 0s
           obj.vals = zeros(size(obj.indxs(obj.indxs(:,1)~=0),1),...
               size(obj.indxs(obj.indxs(:,2)~=0),1),...
               size(obj.indxs(obj.indxs(:,3)~=0),1),...
               size(obj.indxs(obj.indxs(:,4)~=0),1),size(values,2));
           % Loop through values and assign them to the 4 d matrix
           for i=1:size(values,1)
               disp(num2str(i));
               indx1 = obj.indxs(obj.keys(:,1)==data(i,1),1);
               indx1 = indx1(1,1);
               indx2 = obj.indxs(obj.keys(:,2)==data(i,2),2);
               indx2 = indx2(1,1);
               indx3 = obj.indxs(obj.keys(:,3)==data(i,3),3);
               indx3 = indx3(1,1);
               indx4 = obj.indxs(obj.keys(:,4)==data(i,4),4);
               indx4 = indx4(1,1);
               %if there is a value already there, then add it to it, don't
               %overwrite it - this accounts for transshipment
               obj.vals(indx1,indx2,indx3,indx4,:)=obj.vals(indx1,indx2,indx3,indx4,:)...
                   +reshape(values(i,:),1,1,1,1,size(values(i,:),2));
           end %for i
           
           %call function to update for new format
           obj.generateFromHashTable(obj);
        end %function setup
        % ======================================================================
    %> @brief 
    %>
    %> @param 
    %> @param 
    %> @retval
    % ======================================================================
        function setup3d(obj,data,values)
            warning('This will be removed in future versions. Please use the new  format');
           obj.indxs = zeros(size(data,1),size(data,2));
           obj.keys = obj.indxs;
           %Build referencing
           for i=1:size(data,2)
               tmp=sort(unique(data(:,i)),'ascend');
               obj.keys(1:size(tmp,1),i)=tmp;
               obj.indxs(1:size(tmp,1),i) = [1:size(tmp,1)]';
           end %for i
           obj.indxs = obj.indxs(any(obj.indxs,2),:); %remove rows that contain all 0s
           obj.keys = obj.keys(any(obj.keys,2),:); %remove rows that contain all 0s
           obj.vals = zeros(obj.indxs(:,1),obj.indxs(:,2),size(values,2));
           % Loop through values and assign them to the 4 d matrix
           for i=1:size(values,1)
               indx1 = obj.indxs(obj.keys(:,1)==data(i,1),1);
               indx2 = obj.indxs(obj.keys(:,2)==data(i,2),2);
               obj.vals(indx1,indx2,:)=values(i,:);
           end %for i
           %call function to update for new format
           obj.generateFromHashTable(obj);
        end %function setup
        % ======================================================================
    %> @brief 
    %>
    %> @param 
    %> @param 
    %> @retval
    % ======================================================================
        function [x]= getVals5d(obj,indx1,indx2,indx3,indx4,indxsVals)
            warning('This will be removed in future versions. Please use the new  format');
            % indx1 and indx2 and indxsVals should be single vals while indx3 and indx4
            % should be ranges
            tmpSize = size(indx3,1)*size(indx4,1);
            x = reshape(obj.vals(indx1,indx2,indx3,indx4,indxsVals),tmpSize,1);   
        end
         % ======================================================================
    %> @brief 
    %>
    %> @param 
    %> @param 
    %> @retval
    % ======================================================================
        function [x]= getVals3d(obj,indx1,indx2,indxsVals,keys)
            warning('This will be removed in future versions. Please use the new  format');
            if nargin==5
               if keys==2
                   %Loop through the keys and set the indxs
                   keys1 = indx1;
                   keys2 = indx2;
                   indx1=[];
                   indx2=[];
                   %can use ismember as the arrays should be sorted
                   indx1 = ismember(obj.keys(:,1),keys1);
                   indx2 = ismember(obj.keys(:,2),keys2);
               end %if
            end %if
            tmpSize = size(obj.indxs(indx1,1),1)*size(obj.indxs(indx2,2),1);
            tmp = obj.vals;
            tmp(~indx1,:,:)= [];
            tmp(:,~indx2,:)=[];
            tmp = tmp(:,:,indxsVals);
            %x = reshape(obj.vals(indx1,indx2,indxsVals),tmpSize,1);   
            x = reshape(tmp,numel(tmp),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief creates the matrix from matrix data
        %>
        %> @param obj - instance of this class
        %> @param matrix of data coming through. The instance of each row
        %> of keys must be unique
        %> @param valsIndxs array of values corresponding to the index of
        %> dataMatrix that contains keys and values values
        function generateFromTableData(obj,dataMatrix,valsIndxs)
            % keys are stored in the non valsIndxs rows
            obj.keyStore = struct('keys',[]);
            obj.indxStore = struct('indxs',[]);
            obj.noKeys = size(dataMatrix,2) - length(valsIndxs);
            obj.noVals = length(valsIndxs);
            obj.cardinality = zeros(obj.noKeys,1);
            keyIndx = 1;
            %get our keys
            keyIndxs =[];
            disp('Getting the keys');
            for i=1:size(dataMatrix,2)
                if ~ismember(i,valsIndxs)
                    obj.keyStore(keyIndx).keys = full(unique(dataMatrix(:,i)));
                    obj.indxStore(keyIndx).indxs = ...
                        [1:length(obj.keyStore(keyIndx).keys)]';
                    obj.cardinality(keyIndx) =...
                        length(obj.keyStore(keyIndx).keys);
                    
                    keyIndx = keyIndx+1;
                    
                    % Now store the indxs of the key fields
                    keyIndxs = [keyIndxs;i];
                end %if
            end %for i
            
            %now build the sparse data matrix
            noRows = cumprod(obj.cardinality);
            obj.data=sparse(noRows(end),obj.noVals);
            keyVals = dataMatrix(:,keyIndxs);
            keyIndxVals = zeros(size(keyVals));
            
            for j=1:size(keyVals,2)
                [~,indxs] = ismember(keyVals(:,j),obj.keyStore(j).keys);
                keyIndxVals(:,j) = full(indxs);
            end %for j
            indx = GenProb.AssignmentToIndex(keyIndxVals,obj.cardinality);
            for j=1:size(valsIndxs,2)
                disp(sprintf('Filling the data. %d of %d',j, size(valsIndxs,2)));
                obj.data(indx,j) = dataMatrix(:,valsIndxs(j));
            end %for j
%             for i=1:size(dataMatrix,1)
%                 disp(sprintf('Building MultiDDataMapper from table. Row %d of %d',i,size(dataMatrix,1)));
%                 keyVals = dataMatrix(i,keyIndxs);
%                 keyIndxVals = zeros(size(keyVals));
%                 for j=1:length(keyVals)
%                     keyIndxVals(j) = obj.indxStore(j).indxs(...
%                         obj.keyStore(j).keys==keyVals(j));
%                 end %for j
%                 indx = GenProb.AssignmentToIndex(keyIndxVals,obj.cardinality);
%                 for j=1:length(valsIndxs)
%                     obj.data(indx,j) = dataMatrix(i,valsIndxs(j));
%                 end %for j
%             end %for i
        end %generateFromTableData
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Allows backwards compatibility from HashTable to this
        %> function
        %>
        %> @param obj
        %> @param ht instance of hastable that contains the data
        function generateFromHashTable(obj,ht)
            %get the number of dimensions
            obj.noKeys = size(ht.keys,2);
            obj.keys = ht.keys;
            obj.vals = ht.vals;
            obj.keyNames = cell(obj.noKeys,1);
            obj.keyNames = cellfun(@(x) '',obj.keyNames,'un',0);
            obj.noVals = size(ht.vals);
            obj.noVals = obj.noVals(end);
            obj.valNames = cell(obj.noVals,1);
            obj.valNames = cellfun(@(x) '',obj.valNames,'un',0);
            obj.cardinality = zeros(obj.noKeys,1);
            obj.keyStore = struct('keys',[]);
            obj.indxStore = struct('indxs',[]);
            for i=1:obj.noKeys
                obj.cardinality(i) = sum(ht.keys(:,i)~=0);
                obj.keyStore(i).keys = ht.keys(ht.keys(:,i)~=0,i);
                obj.indxStore(i).indxs = [1:length(obj.keyStore(i).keys)]';
            end %for i
            %set the values
            %store the data in table format, not multidimensional
            totalRows = cumprod(obj.cardinality);
            totalRows = totalRows(end);
            obj.data= zeros(totalRows,obj.noVals);
            %Now add the values
            for i=1:obj.noVals
                tmp = ht.vals(:);
                keyCountCard = cumprod(obj.cardinality);
                keyCountCard = keyCountCard(end);
                
                startIndx = keyCountCard*(i-1)+1;
                endIndx = keyCountCard*i;
                obj.data(:,i) = tmp(startIndx:endIndx,:);
                clear('tmp');
            end %for i
            
            %Convert it to a sparse matrix
            obj.data = sparse(obj.data);
            %done
        end %generateFromHashTable
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief outputs the data array but as a dataset and using the
        %> keynames
        %>
        %> @param obj
        %> @param removeZeroRows remove rows where all the vals are 0 
        %> optional [1 is default]
        
        %> @retval ds instance of dataset as described above
        function ds = getDataset(obj,removeZeroRows)
            if nargin==1
                removeZeroRows = 1;
                keysToFilterBy =[];
                keyVals = [];
            end %if
            if nargin == 2 || nargin == 3
                keysToFilterBy =[];
                keyVals = [];
            end %if
            %first make sure keynames actually entered
            outputKeyNames = obj.keyNames;
            for i=1:obj.noKeys
                if strcmp(outputKeyNames{i},'')
                    outputKeyNames{i} = sprintf('key_%d',i);
                end %if
            end %for i
            %repeat for value names
            outputValNames = obj.valNames;
            for i=1:obj.noVals
                if strcmp(outputValNames{i},'')
                    outputValNames{i} = sprintf('val_%d',i);
                end %if
            end %for i
            %now pop in the keys
            ds = dataset;
            keyMatrix = obj.createKeyMatrix();
            for i=1:obj.noKeys
                ds.(outputKeyNames{i}) = obj.keyStore(i).keys(keyMatrix(:,i),:);
            end %for i
            clear('keyMatrix');
            %pop in the vals
            %obj.data = (obj.data);
            for i=1:obj.noVals
                ds.(strrep(strrep(outputValNames{i},' ',''),',','')) = obj.data(:,i);
            end %for i
            %convert the data back to a sparse matrix
            obj.data = sparse(obj.data);
            
            
            if removeZeroRows==1
                %remove rows where the values are all 0
                getZeros = zeros(size(ds,1),obj.noVals);
                for i=1:obj.noVals
                    getZeros(:,i) = double(ds(:,obj.noKeys+i))~=0;
                end %for i
                keep = any(getZeros,2);
                ds = ds(keep,:);
                clear('keep','getZeros');
            end %if
            
        end %function getDataset
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief Gets the indx of the keyname. It assumes that all key
        %> names are unique
        %>
        %> @param obj
        %> @param keyName the name of the key to find (can be a cell array
        %> of strings
        %> @retval the index of the key
        function idx = getIndxOfKey(obj,keyName)
            if isa(keyName,'char')
                idx = cellfun(@(x) strcmpi((x),keyName),obj.keyNames);
                idx = find(idx==1);
            else
                idx = zeros(size(keyName));
                for i=1:length(idx)
                    tmpIdx = cellfun(@(x) strcmpi((x),keyName{i}),obj.keyNames);
                    idx(i) = find(tmpIdx==1);
                end %for i
            end %if
        end %getIdxOfKey
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief Creates a matrix based on the keys - for use with
        %> outputting the dataset
        %>
        %> @param obj
        function keyMatrix = createKeyMatrix(obj)
            totalRows = cumprod(obj.cardinality);
            totalRows = totalRows(end);
            keyMatrix = zeros(totalRows,obj.noKeys);
            for i=1:obj.noKeys
                if i==1
                    noCopiesStack=1;
                    noCopiesSeq = cumprod(obj.cardinality(i+1:end));
                    noCopiesSeq = noCopiesSeq(end);
                elseif i==obj.noKeys
                    noCopiesStack = cumprod(obj.cardinality(1:i-1));
                    noCopiesStack = noCopiesStack(end);
                    noCopiesSeq = 1;
                else
                    noCopiesStack = cumprod(obj.cardinality(1:i-1));
                    noCopiesStack = noCopiesStack(end);
                    noCopiesSeq = cumprod(obj.cardinality(i+1:end));
                    noCopiesSeq = noCopiesSeq(end);
                end %if
                %First copy each value multip times in a row (only for i>1)
                copiesSeq = repmat(obj.indxStore(i).indxs,noCopiesSeq,1);
                %Now stack them (none stacked for last key)
                copiesStack = repmat(copiesSeq',noCopiesStack,1);
                keyMatrix(:,i) = copiesStack(:);
            end %for i
        end %getIdxOfKey
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief gets a dataset of values based on the passed key values
        %>
        %> @param obj
        %> @param returnKeys cell array of key names to return data for
        %> @param returnVals cell array of key names of values to return
        %> @param keysToFilterBy cell array of key names to filter by
        %> @param keyVals structure with an instance for each of
        %> keysToFilterBy. Each structure contains an internal variable
        %> called vals. Each value in vals is an item we want to keep
        %> @retval ds the return dataset
        function ds = getFilteredDataset(obj,returnKeys,returnVals,keysToFilterBy,keyVals,...
                removeZeroRows)
            %get the full dataset first
            ds = obj.getDataset(removeZeroRows);
            %filter by keys
                        
            % Now strip out any data that we don't want to keep
            if ~isempty(keysToFilterBy)
                for i=1:length(keysToFilterBy)
                    ds = ds(ismember(ds.(keysToFilterBy{i}),keyVals(i).vals),:);
                end %for i
            end %if
%             for i=1:length(filterKeys)
%                 if isa(filterKeyVals(i),'double')
%                     ds = ds(ds.(obj.keyNames{filterKeys(i)})==filterKeyVals(i),:);
%                 else
%                     ds = ds(cellfun(@(x) strcmpi(x,filterKeyVals{i}),...
%                         ds.(obj.keyNames{filterKeys(i)})),:);
%                 end %if
%             end %for i
            % Add the value fields that we want to output
            ds_new = ds(:,1:obj.noKeys);
            for i=1:length(returnVals)
                    ds_new.(strrep(strrep(returnVals{i},' ',''),',','')) = ...
                        ds.(strrep(strrep(returnVals{i},' ',''),',',''));
            end %for i
            ds = ds_new;
            clear('ds_new');
            %No integrate over fields that we are removing
            
            %roll through the keys integrating as we go for field we're not
            %keeping
%             returnKeyIndxs = zeros(obj.noKeys,1);
%             for i=1:length(returnKeys)
%                 returnKeyIndxs = returnKeyIndxs + ...
%                     ismember(obj.keyNames,returnKeys(i));
%             end %for 
%             
%            [uniqueRows,dsRows,dsIndxs] = unique(double(ds(:,...
%                logical([returnKeyIndxs;zeros(obj.noVals,1)]))),'rows');
%            %create a new holder for the dataset
%           % valFields = ismember(obj.valNames,returnVals);
%            ds_new = ds(1:length(uniqueRows),...
%                logical([returnKeyIndxs(:);ones(size(returnVals))]));
%            if length(uniqueRows)~=size(ds,1)
%                 %for i=1:length(uniqueRows)
%                     disp(sprintf('Filling dataset%d of %d',i,length(uniqueRows)));
%                     %Fill in the keys for this row
%                     ds_new(i,1:length(returnKeys)) = ds(dsRows(i,:),logical(returnKeyIndxs));
%                     %ds_new(:,1:length(returnKeys)) = ds(dsRows(:,:),logical(returnKeyIndxs));
%                     %now sum over the fields we're aggregating for 
%                     for j=1:length(returnVals)
%                         %j
%                         ds_new.(strrep(strrep(returnVals{j},' ',''),',','')) = ...
%                             sum(double(ds.(strrep(strrep(returnVals{j},' ',''),',',''))(dsIndxs==i,:)));
%                     end %for j
%                 %end %for i
%            end %if
%             ds = ds_new;
            %remove keys we don't want
             if removeZeroRows==1
                %remove rows where the values are all 0
                getZeros = zeros(size(ds,1),length(returnVals));
                for i=1:length(returnVals)
                    getZeros(:,i) = double(ds(:,length(returnKeys)+i))~=0;
                end %for i
                keep = any(getZeros,2);
                ds = ds(keep,:);
                clear('keep','getZeros');
            end %if
            
        end %getFilteredDataset
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief gets the indexes relevant to the passed arrray
        %>
        %> @param obj
        %> @param matrix
        %> @param keys
        %> @param keyVals
        function indxs = getIndxsFromKeys(obj,matrix,keys,keyVals)
            if size(keys(:),2) >1
                indxs = getIndxsFromKeys(obj,matrix,keys(2:end),keyVals(2:end));
            else
                indxs=[];
            end %if
            
            indxs = [matrix(:,keys(1))==keyVals(1),indxs];
        end %getIndxsFromKeys
        %==================================================================
        %> @brief returns the keys given a key name
        %> 
        %> @param obj
        %> @param keyName
        %> 
        %> @retval keyList returns the list of keys
        %> @retval keyIndxs associated indexes with each key
        function [keyList,keyIndxs] = getKeys(obj,keyName)
            keyList = full(obj.keyStore(obj.getIndxOfKey(keyName)).keys);
            keyIndxs = full(obj.indxStore(obj.getIndxOfKey(keyName)).indxs);
        end %function getKeys
    end %methods
    
end

