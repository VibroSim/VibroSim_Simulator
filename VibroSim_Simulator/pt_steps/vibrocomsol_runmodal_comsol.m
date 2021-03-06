% function [ ret ] = runmodal(dc_dest_href,dc_measident_str,dc_model_comsol_href)
import com.comsol.model.util.*

modelfile = dc_model_comsol_href{1};
model=mphload(modelfile);

% Make sure progress bar/plots are displayed
ModelUtil.showProgress(true);
ModelUtil.showPlots(true);

% Run the modal study
model.study('solidmech_modal_study').run;

% Determine filename for modal results

%[modelpath,modelname,modelext] = fileparts(modelfile);
modal_export_name=fullfile(dc_dest_href{1},[ dc_measident_str '_modalfreq.txt' ]);

model.result.export('solidmech_modal_export').set('filename',buildabspath(modal_export_name));
model.result.export('solidmech_modal_export').run();

%% Save .mph output in temporary directory in case use wants to load it back in
%savefilename = fullfile(tempdir,sprintf('vibrosim_%s',char(java.lang.System.getProperty('user.name'))),[ modelname '_modalcalc.mph' ]);

% Save .mph output
savefilename = fullfile(dc_dest_href{1},[ dc_measident_str '_modalcalc.mph' ]);
mphsave(model,savefilename); 

ret={ 
  { 'dc:modalcalc_comsol', {savefilename} },
  { 'dc:modalfreqs', {modal_export_name} },
};


