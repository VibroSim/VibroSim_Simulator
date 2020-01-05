% function [ ret ] = vibrocomsol_multisweep_seg1_analysis_comsol(dc_dest_href,dc_measident_str,dc_model_comsol_withsegboundaries_href,segnum_int)

import com.comsol.model.util.*

modelfile = dc_model_comsol_withsegboundaries_href{1};
model=mphload(modelfile);

% Make sure progress bar/plots are displayed
ModelUtil.showProgress(true);
ModelUtil.showPlots(true);

% Determine filename for harmonic sweep results

%[modelpath,modelname,modelext] = fileparts(modelfile);
%spectrum_image_name=fullfile(modelpath,[ modelname '_sweep_spectrum.png' ]);


% Run the seg harmonic study
model.study(sprintf('solidmech_multisweep_seg%d_study',segnum_int)).run;
% Save .mph output in temporary directory in case use wants to load it back in
savefilename = fullfile(tempdir,sprintf('vibrosim_%s',char(java.lang.System.getProperty('user.name'))),[ dc_measident_str sprintf('_multisweep_seg%d.mph',segnum_int) ]);
mphsave(model,savefilename); 


%% run the plot/export
%model.result('vibro_harmonicsweep_spectrum').run();
%model.result.export('vibro_harmonicsweep_spectrum_image').run()


% Write out the transducer contact probe displacement and velocity
model.result.table(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_displ_table',segnum_int)).clearTableData()
model.result.table(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_displ_table',segnum_int)).set('storetable','inmodelandonfile');

seg_xducercontactprobe_displ_file = fullfile(dc_dest_href{1},sprintf('%s_xducercontactprobe_displspec_seg%d.txt',dc_measident_str,segnum_int));
model.result.table(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_displ_table',segnum_int)).set('filename',buildabspath(seg_xducercontactprobe_displ_file));

model.result.numerical(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_displ',segnum_int)).run;
model.result.numerical(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_displ',segnum_int)).setResult();




model.result.table(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_vel_table',segnum_int)).clearTableData()
model.result.table(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_vel_table',segnum_int)).set('storetable','inmodelandonfile');

seg_xducercontactprobe_vel_file = fullfile(dc_dest_href{1},sprintf('%s_xducercontactprobe_velspec_seg%d.txt',dc_measident_str,segnum_int));
model.result.table(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_vel_table',segnum_int)).set('filename',buildabspath(seg_xducercontactprobe_vel_file));

model.result.numerical(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_vel',segnum_int)).run;
model.result.numerical(sprintf('custom_vibrodynamic_multisweep_seg%d_xducercontactprobe_vel',segnum_int)).setResult();




% Write out the laser displacement and velocity

model.result.table(sprintf('custom_solidmech_multisweep_seg%d_laser_displ_table',segnum_int)).clearTableData()
  model.result.table(sprintf('custom_solidmech_multisweep_seg%d_laser_displ_table',segnum_int)).set('storetable','inmodelandonfile');

seg_laser_displ_file = fullfile(dc_dest_href{1},sprintf('%s_laser_displspec_seg%d.txt',dc_measident_str,segnum_int));
model.result.table(sprintf('custom_solidmech_multisweep_seg%d_laser_displ_table',segnum_int)).set('filename',buildabspath(seg_laser_displ_file));

model.result.numerical(sprintf('custom_solidmech_multisweep_seg%d_laser_displ',segnum_int)).run;
model.result.numerical(sprintf('custom_solidmech_multisweep_seg%d_laser_displ',segnum_int)).setResult();




model.result.table(sprintf('custom_solidmech_multisweep_seg%d_laser_vel_table',segnum_int)).clearTableData()
model.result.table(sprintf('custom_solidmech_multisweep_seg%d_laser_vel_table',segnum_int)).set('storetable','inmodelandonfile');

seg_laser_vel_file = fullfile(dc_dest_href{1},sprintf('%s_laser_velspec_seg%d.txt',dc_measident_str,segnum_int));
model.result.table(sprintf('custom_solidmech_multisweep_seg%d_laser_vel_table',segnum_int)).set('filename',buildabspath(seg_laser_vel_file));

model.result.numerical(sprintf('custom_solidmech_multisweep_seg%d_laser_vel',segnum_int)).run;
model.result.numerical(sprintf('custom_solidmech_multisweep_seg%d_laser_vel',segnum_int)).setResult();


% Write out the crack stress



model.result.table(sprintf('solidmech_multisweep_seg%d_crackcenterstress_table',segnum_int)).clearTableData()
model.result.table(sprintf('solidmech_multisweep_seg%d_crackcenterstress_table',segnum_int)).set('storetable','inmodelandonfile');

seg_crackcenterstress_file = fullfile(dc_dest_href{1},sprintf('%s_crackcenterstressspec_seg%d.txt',dc_measident_str,segnum_int));
model.result.table(sprintf('solidmech_multisweep_seg%d_crackcenterstress_table',segnum_int)).set('filename',buildabspath(seg_crackcenterstress_file));

model.result.numerical(sprintf('solidmech_multisweep_seg%d_crackcenterstress',segnum_int)).run;
model.result.numerical(sprintf('solidmech_multisweep_seg%d_crackcenterstress',segnum_int)).setResult();


ret={ 
  { sprintf('dc:seg%d_xducercontactprobe_displ',segnum_int), {seg_xducercontactprobe_displ_file} },
  { sprintf('dc:seg%d_xducercontactprobe_vel',segnum_int), {seg_xducercontactprobe_vel_file} },
  { sprintf('dc:seg%d_laser_displ',segnum_int), {seg_laser_displ_file} },
  { sprintf('dc:seg%d_laser_vel',segnum_int), {seg_laser_vel_file} },
  { sprintf('dc:seg%d_crackcenterstress',segnum_int), {seg_crackcenterstress_file} },
};

