% function [ ret ] = sweep_analysis(dc_model_comsol_withsweepfrequencies_href)

modelfile = dc_model_comsol_withsweepfrequencies_href{1};
model=mphload(modelfile);

% Make sure progress bar/plots are displayed
ModelUtil.showProgress(true);
ModelUtil.showPlots(true);

% Determine filename for harmonic sweep results

[modelpath,modelname,modelext] = fileparts(modelfile);
spectrum_image_name=fullfile(modelpath,[ modelname '_sweep_spectrum.png' ]);


CreateOrReplace(model.result.export,'vibro_harmonicsweep_spectrum_image','vibro_harmonicsweep_spectrum','Image1D');
model.result.export('vibro_harmonicsweep_spectrum_image').label('vibro_harmonicsweep_spectrum_image');
model.result.export('vibro_harmonicsweep_spectrum_image').set('pngfilename',buildabspath(spectrum_image_name));
%model.result.export('vibro_harmonicsweep_spectrum_image').run();


% Run the harmonic study
model.study('solidmech_harmonicsweep_study').run;


% Save .mph output in temporary directory in case use wants to load it back in
savefilename = fullfile(tempdir,sprintf('vibrosim_%s',char(java.lang.System.getProperty('user.name'))),[ modelname '_harmonicsweep.mph' ]);
mphsave(model,savefilename); 

% run the plot/export
model.result('vibro_harmonicsweep_spectrum').run();
model.result.export('vibro_harmonicsweep_spectrum_image').run();

ret={ 
  { 'dc:sweep_spectrum', {spectrum_image_name} },
};

