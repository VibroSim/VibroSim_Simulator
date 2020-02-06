% function [ ret ] = setrayleighfactors(dc_dest_href,dc_measident_str,dc_model_comsol_href,dc_spcrayleighdamping_alpha_float,dc_spcrayleighdamping_beta_float)

import com.comsol.model.util.*

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

model.param.set('spcrayleighdamping_alpha',[ to_string(dc_spcrayleighdamping_alpha_float) '[s^-1]' ]);
model.param.set('spcrayleighdamping_beta',[ to_string(dc_spcrayleighdamping_beta_float) '[s]' ]);

%[modelpath,modelname,modelext] = fileparts(modelfile);

% model_comsol_withsweepfrequencies_name = fullfile(dc_dest_href{1},[ dc_measident_str '_withrayleigh.mph' ]);


mphsave(model,modelfile); 


ret={ 
  { 'dc:model_comsol', {modelfile} }
};

