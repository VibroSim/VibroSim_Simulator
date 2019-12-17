% function [ ret ] = heatflow_analysis(dc_dest_href,dc_measident_str,dc_model_comsol_href,dc_heatingdata_href,dc_exc_t3_float)

import com.comsol.model.util.*

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

% Make sure progress bar/plots are displayed
ModelUtil.showProgress(true);
ModelUtil.showPlots(true);

% refresh the crack heating function data (from the previous step) from disk
%model.func('crack_heatingfunction').refresh;
model.func('crack_heatingfunction').set('filename',buildabspath(dc_heatingdata_href{1})); % must use absolute path here because COMSOL current directory is the mli (LiveLink) directory


% Determine filename for heating results

%[modelpath,modelname,modelext] = fileparts(modelfile);
heating_image_name=fullfile(dc_dest_href{1},[ dc_measident_str '_heating.png' ]);

CreateOrReplace(model.result.export,'vibro_heating_image','vibro_heating_plot','Image3D');
model.result.export('vibro_heating_image').label('vibro_heating_image');
model.result.export('vibro_heating_image').set('pngfilename',buildabspath(heating_image_name));
%model.result.export('vibro_heating_image').run();


% Run the heatflow study
model.study('heatflow_study').run;


% update the solution index (which time) for the plot
vibro_heating_dset=model.result('vibro_heating_plot').getString('data');
vibro_heating_solution=model.result.dataset(vibro_heating_dset).getString('solution');
vibro_heating_times=model.sol(vibro_heating_solution).getPVals();

% find the time at or immediately prior to t3
times_rel_t3 = vibro_heating_times-dc_exc_t3_float;
timeindexes_through_t3 = find(times_rel_t3 <= 0.0);
plot_index = timeindexes_through_t3(length(timeindexes_through_t3));

model.result('vibro_heating_plot').set('solnum',plot_index);
% Run the plot/export
model.result('vibro_heating_plot').run();
model.result.export('vibro_heating_image').run();

% Save .mph output
savefilename = fullfile(dc_dest_href{1},[ dc_measident_str '_heating.mph' ]);
mphsave(model,savefilename); 



ret={ 
  { 'dc:vibro_heating_image', {heating_image_name} },
  { 'dc:heatflow_comsol', {savefilename} },
};
