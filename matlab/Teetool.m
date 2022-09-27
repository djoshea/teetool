classdef Teetool < handle
    % provides matlab interface to Teetool by Willem Eerland:
    % paper: https://doi.org/10.5334/jors.163
    % code: https://doi.org/10.5281/zenodo.251481, 
    % using Python3 version at https://github.com/sfo/teetool
    
    methods(Static)
        function world = build_world(timeCell, dataCell, args)
            % dataCell is nConditions { nTrials x time x d } where d is 2 or 3
            arguments
                timeCell
                dataCell cell
                args.resolution (1, :) = [100 100];
                args.model_type (1, 1) string = "resampling";
                args.ngaus (1, 1) = 100;
            end
        
            nC = numel(dataCell);
            D = size(dataCell{1}, 3);

            assert(D == 2 || D == 3, "D==%d is not valid", D);
            world = py.teetool.World(name="world", ndim=int32(D), resolution=py.list(uint32(args.resolution)));

            prog = ProgressBar(nC, "Adding conditions to teetool.World");
            %teetool_data = cell(nC, 1);
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
        
                this_cond = cell(1, R);
                for iR = 1:R
                    % create list of tuples (time, data) for each trial
                    this_time = shiftdim(time(iR, :), 1); % T x 1
                    this_data = shiftdim(data(iR, :, :), 1); % T x D
                    tmask = ~isnan(this_time) & ~any(isnan(this_data), 2);
                    this_cond{iR} = py.tuple({this_time(tmask), this_data(tmask, :)});
                end
        
                this_input = py.list(this_cond);
                world.addCluster(this_input, "c" + iC);
                prog.update(iC);
            end
            prog.finish();
            
            debug('Building world model...');

            settings = py.dict(model_type=args.model_type, ngaus=int32(args.ngaus));
            world.buildModel(settings);
            fprintf('done\n');
        end

        function out = get_log_likelihood_data(world)
            arguments
                world py.teetool.world.World
            end
            % out is nC x 1 struct, for each condition:
            %   .xvec and .yvec are the x and y values to pass to meshgrid
            %   .logp is the loglikelihood values
            nC = length(py.getattr(world, '_clusters'));

            prog = ProgressBar(nC, 'Extracting log likelihood images');
            for iC = nC:-1:1
                logp = world.getLogLikelihood(py.list({int32(iC-1)}));
                ll = single(logp{1}{1})';
                xx = single(logp{2}{1})';
                yy = single(logp{2}{2})'; 
                out(iC).logp = ll;
                out(iC).xvec = xx(1, :);
                out(iC).yvec = yy(:, 1);
                prog.increment();
            end
            prog.finish();
            out = out';
        end

        function out = get_sd_contours(world, sd_widths)
            arguments
                world py.teetool.world.World
                sd_widths (:, 1) = [1; norminv(1-0.05/2); norminv(1-0.01/2)];
            end
            % out is nC x numel(sdwidths) struct, for each condition:
            %   .xvec and .yvec are the x and y values to pass to meshgrid
            %   .logp is the loglikelihood values
            nC = length(py.getattr(world, '_clusters'));
            nS = numel(sd_widths);

            prog = ProgressBar(nC * nS, 'Extracting contour hulls');
            figh = figure(Visible=false);
            axh = axes('Parent', figh);
            for iC = nC:-1:1
                for iS = 1:nS
                    tube = world.getTube(py.list({int32(iC-1)}), sdwidth=sd_widths(iS));
                    inside = logical(tube{1}{1})';
                    xx = double(tube{2}{1})';
                    yy = double(tube{2}{2})';
                    
                    [cv, h] = contour(xx, yy, inside, [1 1], Parent=axh);
                    delete(h);
                    
                    out(iC, iS).inside = sparse(inside);
                    [out(iC, iS).xvecs, out(iC, iS).yvecs] = parse_countours(cv);
                    if numel(out(iC, iS).xvecs) > 1
                        a = 1;
                    end
                    out(iC, iS).sd_width = sd_widths(iS);
                    prog.increment();
                end
            end

            close(figh);
            prog.finish();

            function [xc, yc] = parse_countours(cv)
                % parse successive contours
                cstart = 1;
                ctotal = size(cv, 2);
                piece_idx = 1;
                xc = {};
                yc = {};
                while cstart <= ctotal
                    npts = cv(2, cstart);
                    xc{piece_idx} = cv(1, cstart + (1:npts)); %#ok<AGROW> 
                    yc{piece_idx} = cv(2, cstart + (1:npts)); %#ok<AGROW> 
                    cstart = cstart + npts + 1;
                    piece_idx = piece_idx + 1;
                end
            end
        end

        function out = get_means(world)
            arguments
                world py.teetool.world.World
            end
            % out is nC x 1 struct, for each condition:
            %   .xvec and .yvec are the x and y values to pass to meshgrid
            nC = length(py.getattr(world, '_clusters'));

            prog = ProgressBar(nC, 'Extracting mean trajectories');
            for iC = nC:-1:1
                m = world.getMean(py.list({int32(iC-1)}));
                xyz = single(m{1});
                out(iC).xvec = xyz(:, 1);
                out(iC).yvec = xyz(:, 2);
                if size(xyz, 2) > 2
                    out(iC).zvec = xyz(:, 3);
                end
                prog.increment();
            end
            prog.finish();
            out = out';
        end
    end
end
