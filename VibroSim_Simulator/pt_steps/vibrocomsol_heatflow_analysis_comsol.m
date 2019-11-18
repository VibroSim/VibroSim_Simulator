% function [ ret ] = heatflow_analysis(dc_model_comsol_href)

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

% Make sure progress bar/plots are displayed
ModelUtil.showProgress(true);
ModelUtil.showPlots(true);

% refresh the crack heating function data (from the previous step) from disk
model.func('crack_heatingfunction').refresh;


% Determine filename for heating results

[modelpath,modelname,modelext] = fileparts(modelfile);
heating_image_name=fullfile(modelpath,[ modelname '_heating.png' ]);

CreateOrReplace(model.result.export,'vibro_heating_image','vibro_heating_plot','Image3D');
model.result.export('vibro_heating_image').label('vibro_heating_image');
model.result.export('vibro_heating_image').set('pngfilename',buildabspath(heating_image_name));
%model.result.export('vibro_heating_image').run();


% Run the heatflow study
model.study('heatflow_study').run;

% Save .mph output
savefilename = fullfile(modelpath,[ modelname '_heating.mph' ]);
mphsave(model,savefilename); 

% Run the plot/export
model.result('vibro_heating_plot').run();
model.result.export('vibro_heating_image').run();


ret={ 
  { 'dc:vibro_heating_image', {heating_image_name} },
  { 'dc:heatflow_comsol', {savefilename} },
};
