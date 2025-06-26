classdef JuliaData < handle
    properties (Access = public)
        nFeatures
        nSamples
        nLabels
        Xtrain
        Ytrain
        Xtest
        Ytest
        Ntest
        batchSize
        Batch_nD
        Batch_nB
        muX
        sigmaX
        muY
        sigmaY
    end

    properties (Access = private)
        data  % Struct returned by Julia constructor with all properties
    end

    methods (Access = public)

        function obj = JuliaData(params)
            % Call Julia constructor
            obj.data = callJuliaClass('Data', 'DataStruct', params);

            % Assign public properties
            obj.nFeatures = obj.data.nFeatures;
            obj.nSamples  = obj.data.nSamples;
            obj.nLabels   = obj.data.nLabels;
            size(obj.data.Xtrain)

            obj.Xtrain = obj.data.Xtrain;
            obj.Ytrain = obj.data.Ytrain;
            obj.Xtest = obj.data.Xtest;
            obj.Ytest = obj.data.Ytest;
            %{
            obj.Xtrain    = cell2mat(obj.data.Xtrain);
            obj.Ytrain    = cell2mat(obj.data.Ytrain);
            obj.Xtest     = cell2mat(obj.data.Xtest);
            obj.Ytest     = cell2mat(obj.data.Ytest);
            %}
            obj.Ntest     = obj.data.Ntest;
            obj.batchSize = obj.data.batchSize;
            obj.Batch_nD  = obj.data.Batch_nD;
            obj.Batch_nB  = obj.data.Batch_nB;
            obj.muX       = obj.data.muX(:)'; % ATTENTION
            obj.sigmaX    = obj.data.sigmaX(:)'; % ATTENTION
            obj.muY       = obj.data.muY(:)'; % ATTENTION
            obj.sigmaY    = obj.data.sigmaY(:)'; % ATTENTION
        end

    end


end