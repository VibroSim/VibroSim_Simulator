% function [ ret ] = setburstfrequency(dc_model_comsol_href,dc_excitation_frequency_float)

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

model.param.set('simulationburstfreq',[ to_string(dc_excitation_frequency_float) '[Hz]' ]);

[modelpath,modelname,modelext] = fileparts(modelfile);

model_comsol_withburstfrequency_name = fullfile(modelpath,[ modelname '_withburstfrequency.mph' ]);


mphsave(model,model_comsol_withburstfrequency_name); 


ret={ 
  { 'dc:model_comsol_withburstfrequency', {model_comsol_withburstfrequency_name} }
};


