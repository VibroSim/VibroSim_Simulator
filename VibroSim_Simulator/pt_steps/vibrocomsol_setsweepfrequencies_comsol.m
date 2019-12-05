% function [ ret ] = setsweepfrequencies(dc_dest_href,dc_measident_str,dc_model_comsol_href,dc_sweep_start_frequency_float,dc_sweep_step_frequency_float,dc_sweep_end_frequency_float)

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

model.param.set('simulationfreqstart',[ to_string(dc_sweep_start_frequency_float) '[Hz]' ]);
model.param.set('simulationfreqstep',[ to_string(dc_sweep_step_frequency_float) '[Hz]' ]);
model.param.set('simulationfreqend',[ to_string(dc_sweep_end_frequency_float) '[Hz]' ]);

%[modelpath,modelname,modelext] = fileparts(modelfile);

model_comsol_withsweepfrequencies_name = fullfile(dc_dest_href{1},[ dc_measident_str '_withsweepfrequencies.mph' ]);


mphsave(model,model_comsol_withsweepfrequencies_name); 


ret={ 
  { 'dc:model_comsol_withsweepfrequencies', {model_comsol_withsweepfrequencies_name} }
};


