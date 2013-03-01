%> @file HashTable.m
%> @brief Simple class for creating a hashtable
%> @section matlabComments Details
%> @authors Eoin O'Keeffe (eoin.okeeffe.09@ucl.ac.uk)
%> @date initiated: 17/08/2011
%> @version 
%> 1.0
%> @section intro Method
%> Takes in a multivariable array, and indexes it so it can be referenced
%> properly. Created for IEA for the large datasets
%> @subsection version_history
%> 1.0 Initial set up for IEA
%>@attention
%>@todo 
classdef HashTable < handle
    %HASHTABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        indxs
        keys
        vals
    end
    
    methods
        % ======================================================================
    %> @brief 
    %>
    %> @param 
    %> @param 
    %> @retval
    % ======================================================================
        function setup5d(obj,data,values)
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
        end %function setup
        % ======================================================================
    %> @brief 
    %>
    %> @param 
    %> @param 
    %> @retval
    % ======================================================================
        function setup3d(obj,data,values)
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
        end %function setup
        % ======================================================================
    %> @brief 
    %>
    %> @param 
    %> @param 
    %> @retval
    % ======================================================================
        function [x]= getVals5d(obj,indx1,indx2,indx3,indx4,indxsVals)
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
            x = reshape(obj.vals(indx1,indx2,indxsVals),tmpSize,1);   
        end
        
    end %methods
    
end

