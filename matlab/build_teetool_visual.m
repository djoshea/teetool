classdef Teetool < handle
    methods(Static)
        function teetool_data = build_teetool_visual(timeCell, dataCell, args)
        % dataCell is nConditions { nTrials x d } where d is 2 or 3
        arguments
            timeCell
            dataCell cell
        end
    
        nC = numel(dataCell);
        D = NaN;
        teetool_data = cell(nC, 1);
        for iC = 1:nC
            data = dataCell{iC};
    
            [R, T, thisD] = size(data);
            
            if iscell(timeCell)
                time = timeCell{iC};
                assert(size(time, 1) == R);
                assert(size(time, 2) == T);
            else
                assert(numel(timeCell) == T);
                time = repmat(makerow(timeCell), R, 1);
            end
    
            if isnan(D)
                D = thisD;
            else
                assert(thisD == D);
            end
    
            this_cond = cellvec(R);
            for iR = 1:R
                % create list of tuples (time, data) for each trial
                this_cond{iR} = py.tuple(shiftdim(time(iR, :), 1), shiftdim(data(iR, :, :), 1));
            end
    
            teetool_data{iC} = py.list(this_cond);
        end
        end
end
