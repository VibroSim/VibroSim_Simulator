% function [ ret ] = burst_analysis(dc_dest_href,dc_measident_str,dc_model_comsol_withburstfrequency_href)

modelfile = dc_model_comsol_withburstfrequency_href{1};
model=mphload(modelfile);

% Make sure progress bar/plots are displayed
ModelUtil.showProgress(true);
ModelUtil.showPlots(true);

% Run the harmonic burst study
model.study('solidmech_harmonicburst_study').run;

[normalstrain,shearstrain]=CrackStrain(model,'Geom','crack','solidmech_harmonicburst',1);

% Determine filename for harmonic results

%[modelpath,modelname,modelext] = fileparts(modelfile);
%modal_export_name=fullfile(modelpath,[ modelname '_modalfreq.txt' ]);

%model.result.export('solidmech_modal_export').set('filename',buildabspath(modal_export_name));
%model.result.export('solidmech_modal_export').run();

%% Save .mph output in temporary directory in case use wants to load it back in
%savefilename = fullfile(tempdir,sprintf('vibrosim_%s',char(java.lang.System.getProperty('user.name'))),[ modelname '_burstcalc.mph' ]);

% Save .mph output
savefilename = fullfile(dc_dest_href{1},[ dc_measident_str '_burstcalc.mph' ]);

mphsave(model,savefilename); 

ret={ 
  { 'dc:burstcalc_comsol', {savefilename} },
  { 'dc:harmonicburst_normalstrain', normalstrain },
  { 'dc:harmonicburst_shearstrain', shearstrain },
};

