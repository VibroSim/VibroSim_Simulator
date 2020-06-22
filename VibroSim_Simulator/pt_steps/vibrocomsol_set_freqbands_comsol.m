% function [ ret ] = vibrocomsol_set_freqbands_comsol(dc_dest_href,dc_measident_str,dc_model_comsol_href,dc_freqband_seg1_start_float,dc_freqband_seg1_step_float,dc_freqband_seg1_end_float,dc_freqband_seg2_start_float,dc_freqband_seg2_step_float,dc_freqband_seg2_end_float,dc_freqband_seg3_start_float,dc_freqband_seg3_step_float,dc_freqband_seg3_end_float,dc_freqband_seg4_start_float,dc_freqband_seg4_step_float,dc_freqband_seg4_end_float)

import com.comsol.model.util.*

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

model.param.set('seg1_freqstart',[ to_string(dc_freqband_seg1_start_float) '[Hz]']);
model.param.set('seg1_freqstep',[ to_string(dc_freqband_seg1_step_float) '[Hz]']);
model.param.set('seg1_freqend',[ to_string(dc_freqband_seg1_end_float) '[Hz]']);

model.param.set('seg2_freqstart',[ to_string(dc_freqband_seg2_start_float) '[Hz]']);
model.param.set('seg2_freqstep',[ to_string(dc_freqband_seg2_step_float) '[Hz]']);
model.param.set('seg2_freqend',[ to_string(dc_freqband_seg2_end_float) '[Hz]']);

model.param.set('seg3_freqstart',[ to_string(dc_freqband_seg3_start_float) '[Hz]']);
model.param.set('seg3_freqstep',[ to_string(dc_freqband_seg3_step_float) '[Hz]']);
model.param.set('seg3_freqend',[ to_string(dc_freqband_seg3_end_float) '[Hz]']);

model.param.set('seg4_freqstart',[ to_string(dc_freqband_seg4_start_float) '[Hz]']);
model.param.set('seg4_freqstep',[ to_string(dc_freqband_seg4_step_float) '[Hz]']);
model.param.set('seg4_freqend',[ to_string(dc_freqband_seg4_end_float) '[Hz]']);


%[modelpath,modelname,modelext] = fileparts(modelfile);

model_comsol_withsegboundaries_name = fullfile(dc_dest_href{1},[ dc_measident_str '_withsegboundaries.mph' ]);


mphsave(model,model_comsol_withsegboundaries_name); 


ret={ 
  { 'dc:model_comsol_withsegboundaries', {model_comsol_withsegboundaries_name} }
};


