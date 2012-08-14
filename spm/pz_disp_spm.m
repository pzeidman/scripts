function pz_disp_spm(SPM)
    % Displays trials and contrasts of an SPM GLM for inspection.
    %
    % To use, load your SPM then type: dispSPM(SPM)
    %
    % ---------------------------------------------------------------------
    % Copyright (C) 2012 Peter Zeidman
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses/>.   
    % ---------------------------------------------------------------------        
         
    display_trials();
    display_contrasts();
    
    % ---------------------------------------------------------------------            
    function display_trials()
        
        % We buffer output to show sessions side-by-side
        output = {};
           
        for sess = 1:length(SPM.Sess)
            output{1,sess} = sprintf('SESSION %d',sess);

            % List all trials for this session
            all_onsets = [];
            all_durations = [];
            all_names  = {};
            for i = 1:length(SPM.Sess(sess).U)
                num_trials = length(SPM.Sess(sess).U(i).ons);

                % Get the column name, ignore parametrics
                column_name = SPM.Sess(sess).U(i).name;
                if length(column_name) > 1
                    column_name = column_name{1};
                end

                % Ensure onsets and durations are column vectors
                onset_vector = SPM.Sess(sess).U(i).ons;
                duration_vector = SPM.Sess(sess).U(i).dur;
                if isrow(onset_vector)
                    onset_vector = onset_vector';
                end
                if isrow(duration_vector)
                    duration_vector = duration_vector';
                end
                
                % Record
                all_onsets = vertcat(all_onsets, onset_vector);
                all_names  = vertcat(all_names,  cellstr(repmat(column_name,num_trials,1)));
                all_durations = vertcat(all_durations, duration_vector);
            end

            % Sort names & durations by onsets
            [sorted_onsets, sorted_onsets_idx] = sort(all_onsets);
            sorted_names     = all_names(sorted_onsets_idx);
            sorted_durations = all_durations(sorted_onsets_idx);

            % Build table for session
            output{2,sess} = sprintf('=================================================');
            output{3,sess} = sprintf('%5s %20s %10s %10s','Trial','Column','Onset','Duration');
            output{4,sess} = sprintf('=================================================');
            for i = 1:length(sorted_names)

                % Marker for this trial
                if i > 1 && (sorted_onsets(i) == sorted_onsets(i-1))
                    % Same onset as previous
                    marker = '[';
                elseif i < length(sorted_names) && (sorted_onsets(i) == sorted_onsets(i+1))
                    % Same onset as next
                    marker = '[';
                else
                    % Misc
                    marker = '';
                end

                output{5+i-1,sess} = sprintf('%2s %5d %20s %8.2f %8.2f',marker, i, sorted_names{i}, sorted_onsets(i), sorted_durations(i));
            end        
        end % End session

        % Display
        for row = 1:size(output,1)
            for col = 1:size(output,2)
                fprintf('%50s \t',output{row,col});
            end
            fprintf('\n');
        end
    end
    
    % ---------------------------------------------------------------------
    function display_contrasts()
        fprintf('=================================================\n');
        fprintf('Contrasts\n');
        fprintf('=================================================\n');   

        if ~isfield(SPM,'xCon')
            fprintf('None\n');
            return;
        end
        
        for i = 1:length(SPM.xCon)
           cprintf('-text','%d. %s\n', i, SPM.xCon(i).name);

           % Get non-zero columns from contrast
           nzcolumns = find(SPM.xCon(i).c);
           if strcmp(SPM.xCon(i).STAT,'T')
                for j = 1:length(nzcolumns)
                    contrast_entry = SPM.xCon(i).c(nzcolumns(j));
                    % Choose colour
                    if contrast_entry > 0
                        colour = [0 .5 0];
                    elseif contrast_entry < 0
                        colour = [0.5 0 0];
                    else
                        colour = 'black';
                    end
                    % Colour - not compatible with Matlab 2012
                    % cprintf(colour, '\t %2d %s \n', contrast_entry, SPM.xX.name{nzcolumns(j)});
                    % BW
                    fprintf('\t %2d %s \n', contrast_entry, SPM.xX.name{nzcolumns(j)});
                end
           end

       end
    end
end