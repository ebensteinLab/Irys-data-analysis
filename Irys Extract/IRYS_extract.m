function varargout = IRYS_extract(varargin)
% IRYS_EXTRACT MATLAB code for IRYS_extract.fig
%      IRYS_EXTRACT, by itself, creates a new IRYS_EXTRACT or raises the existing
%      singleton*.
%
%      H = IRYS_EXTRACT returns the handle to a new IRYS_EXTRACT or the handle to
%      the existing singleton*.
%
%      IRYS_EXTRACT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IRYS_EXTRACT.M with the given input arguments.
%
%      IRYS_EXTRACT('Property','Value',...) creates a new IRYS_EXTRACT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IRYS_extract_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IRYS_extract_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IRYS_extract

% Last Modified by GUIDE v2.5 10-Jan-2017 14:53:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IRYS_extract_OpeningFcn, ...
                   'gui_OutputFcn',  @IRYS_extract_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before IRYS_extract is made visible.
function IRYS_extract_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IRYS_extract (see VARARGIN)

% Choose default command line output for IRYS_extract
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IRYS_extract wait for user response (see UIRESUME)
% uiwait(handles.figure1);
curr_dir=pwd;
if ~isempty(dir([curr_dir,'\parameters.mat']))
    set(handles.settings_path,'String',[curr_dir,'\parameters.mat']);
    load_settings_Callback(hObject, eventdata, handles);
end

if ~isempty(dir([curr_dir,'\interest_zones.mat']))
    load([curr_dir,'\interest_zones.mat']);
    interest_zones2=cell(size(interest_zones,1),1);
    for i1=1:size(interest_zones,1)
        zone_to_add=[num2str(interest_zones(i1,1),'%d'),',',num2str(interest_zones(i1,2),'%d'),',',num2str(interest_zones(i1,3),'%d')];
        interest_zones2{i1}=zone_to_add;
    end
    set(handles.interest_zones,'String',interest_zones2);
end

% --- Outputs from this function are returned to the command line.
function varargout = IRYS_extract_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function directory_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function directory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmap_directory_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function xmap_directory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in scans_directories.
function scans_directories_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function scans_directories_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function confidence_lower_limit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function confidence_lower_limit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extra_length_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function extra_length_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xedges_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function xedges_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sum_edges_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sum_edges_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
global i3 mapping_limits mapping_summary averaged_molecules;

try
    settings.directory = get(handles.directory, 'String');
    settings.sum_edges=str2double(get(handles.sum_edges,'String'));
    settings.xedges=str2double(get(handles.xedges,'String'));
    settings.extra_length=str2double(get(handles.extra_length,'String'));
    settings.confidence_lower_limit=str2double(get(handles.confidence_lower_limit,'String'));
    settings.xmap_directory = get(handles.xmap_directory, 'String');
    settings.build_database=false;%logical(get(handles.build_database,'Value'));
    settings.build_mol_list=false;%logical(get(handles.build_mol_list,'Value'));
    settings.scans_directories=get(handles.scans_directories,'String');
    settings.bed_res=str2double(get(handles.bed_res,'String'));
    settings.continue_past=logical(get(handles.continue_past,'Value'));
    settings.write_bed=logical(get(handles.write_bed,'Value'));
    settings.write_graph=logical(get(handles.write_graph,'Value'));
    settings.write_images=logical(get(handles.write_images,'Value'));
    settings.write_images_FOV=logical(get(handles.writeFOVs,'Value'));
    settings.chr_lengths=get(handles.chr_lengths,'data');
    if size(settings.chr_lengths,2)~=1
        settings.chr_lengths=cell2mat(settings.chr_lengths(:,1));
    end
    settings.background_samples_per_pixel=str2double(get(handles.background_samples_per_pixel,'String'));
    settings.background_sub_max_images=str2double(get(handles.background_sub_max_images,'String'));
    settings.peak_detection_threshold=str2double(get(handles.peak_detection_threshold,'String'));
    settings.min_num_of_peaks=str2double(get(handles.min_num_of_peaks,'String'));
    settings.good_fit_cond=str2double(get(handles.good_fit_cond,'String'));
    settings.mapping_range_extra=str2double(get(handles.mapping_range_extra,'String'));
    settings.autosave_time=str2double(get(handles.autosave_time,'String'));
    settings.hmc_threshold=str2double(get(handles.hmc_threshold,'String'));
    settings.hmc_num_peaks=str2double(get(handles.hmc_num_peaks,'String'));
    settings.min_mol_length=str2double(get(handles.min_mol_length,'String'));
    settings.interest_zones=get(handles.interest_zones,'String');
    settings.interest_zones_check=logical(get(handles.interest_zones_check,'Value'));
    settings.selected_IDs_check=logical(get(handles.selected_IDs_check,'Value'));
    settings.selected_IDs_list=get(handles.selected_IDs_list,'String');
    settings.align_ch_num=logical(get(handles.align_ch_num,'Value'));
    settings.include_non_mapped=logical(get(handles.include_non_mapped,'Value'));
    if settings.directory(end)~='\'
        settings.directory=[settings.directory,'\'];
    end
    if settings.xmap_directory(end)~='\'
        settings.xmap_directory=[settings.xmap_directory,'\'];
    end
    warning('off','MATLAB:iofun:UnsupportedEncoding');
    cla(handles.axes1);
    set(handles.status_bar,'String','Starting...');
    drawnow;
    xmap_irys_method_func(settings,handles);
catch
    set(handles.status_bar,'String',['Error']);
end

% --- Executes on button press in continue_past.
function continue_past_Callback(hObject, eventdata, handles)

% --- Executes on button press in write_bed.
function write_bed_Callback(hObject, eventdata, handles)

% --- Executes on button press in write_graph.
function write_graph_Callback(hObject, eventdata, handles)

% --- Executes on button press in write_images.
function write_images_Callback(hObject, eventdata, handles)

% --- Executes on button press in interest_zones_check.
function interest_zones_check_Callback(hObject, eventdata, handles)

function interest_chr_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function interest_chr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interest_start_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function interest_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interest_end_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function interest_end_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bed_res_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function bed_res_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save_stop.
function save_stop_Callback(hObject, eventdata, handles)

% --- Executes when entered data in editable cell(s) in chr_lengths.
function chr_lengths_CellEditCallback(hObject, eventdata, handles)

function scan_dir_to_add_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function scan_dir_to_add_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in add_scan_dir.
function add_scan_dir_Callback(hObject, eventdata, handles)
scans_directories=get(handles.scans_directories,'String');
dir_to_add=get(handles.scan_dir_to_add,'String');
if dir_to_add(end)~='\'
    dir_to_add=[dir_to_add,'\'];
end
scans_directories{numel(scans_directories)+1}=dir_to_add;
set(handles.scans_directories,'String',scans_directories);

% --- Executes on button press in delete_scan_dir.
function delete_scan_dir_Callback(hObject, eventdata, handles)
scans_directories=get(handles.scans_directories,'String');
if length(scans_directories)>1
    scans_directories(get(handles.scans_directories,'value'))=[];
    set(handles.scans_directories,'String',scans_directories,'value',length(scans_directories));
end

% --- Executes on button press in load_settings.
function load_settings_Callback(hObject, eventdata, handles)
load(get(handles.settings_path, 'String'));
set(handles.directory,'String',settings.directory);
set(handles.sum_edges,'String',num2str(settings.sum_edges,'%d'));
set(handles.xedges,'String',num2str(settings.xedges,'%d'));
set(handles.extra_length,'String',num2str(settings.extra_length,'%.2f'));
set(handles.confidence_lower_limit,'String',num2str(settings.confidence_lower_limit,'%.2f'));
set(handles.xmap_directory,'String',settings.xmap_directory);
set(handles.scans_directories,'String',settings.scans_directories);
set(handles.bed_res,'String',num2str(settings.bed_res,'%d'));
set(handles.continue_past,'Value',settings.continue_past);
set(handles.write_bed,'Value',settings.write_bed);
set(handles.write_graph,'Value',settings.write_graph);
set(handles.write_images,'Value',settings.write_images);
set(handles.writeFOVs,'Value',settings.write_images_FOV);
set(handles.chr_lengths,'data',settings.chr_lengths);
set(handles.background_samples_per_pixel,'String',num2str(settings.background_samples_per_pixel,'%d'));
set(handles.background_sub_max_images,'String',num2str(settings.background_sub_max_images,'%d'));
set(handles.peak_detection_threshold,'String',num2str(settings.peak_detection_threshold,'%d'));
set(handles.min_num_of_peaks,'String',num2str(settings.min_num_of_peaks,'%d'));
set(handles.good_fit_cond,'String',num2str(settings.good_fit_cond,'%d'));
set(handles.mapping_range_extra,'String',num2str(settings.mapping_range_extra,'%d'));
set(handles.autosave_time,'String',num2str(settings.autosave_time,'%d'));
set(handles.hmc_threshold,'String',num2str(settings.hmc_threshold,'%d'));
set(handles.hmc_num_peaks,'String',num2str(settings.hmc_num_peaks,'%d'));
set(handles.min_mol_length,'String',num2str(settings.min_mol_length,'%d'));
set(handles.interest_zones_check,'Value',settings.interest_zones_check);
set(handles.interest_zones,'String',settings.interest_zones);
set(handles.selected_IDs_check,'Value',settings.selected_IDs_check);
set(handles.selected_IDs_list,'String',settings.selected_IDs_list);
set(handles.align_ch_num,'Value',settings.align_ch_num);
set(handles.include_non_mapped,'Value',settings.include_non_mapped);

function settings_path_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function settings_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save_settings.
function save_settings_Callback(hObject, eventdata, handles)
settings.directory = get(handles.directory, 'String');
settings.sum_edges=str2double(get(handles.sum_edges,'String'));
settings.xedges=str2double(get(handles.xedges,'String'));
settings.extra_length=str2double(get(handles.extra_length,'String'));
settings.confidence_lower_limit=str2double(get(handles.confidence_lower_limit,'String'));
settings.xmap_directory = get(handles.xmap_directory, 'String');
settings.scans_directories=get(handles.scans_directories,'String');
settings.bed_res=str2double(get(handles.bed_res,'String'));
settings.continue_past=logical(get(handles.continue_past,'Value'));
settings.write_bed=logical(get(handles.write_bed,'Value'));
settings.write_graph=logical(get(handles.write_graph,'Value'));
settings.write_images=logical(get(handles.write_images,'Value'));
settings.write_images_FOV=logical(get(handles.writeFOVs,'Value'));
settings.chr_lengths=get(handles.chr_lengths,'data');
if size(settings.chr_lengths,2)~=1
    settings.chr_lengths=cell2mat(settings.chr_lengths(:,1));
end
settings.background_samples_per_pixel=str2double(get(handles.background_samples_per_pixel,'String'));
settings.background_sub_max_images=str2double(get(handles.background_sub_max_images,'String'));
settings.peak_detection_threshold=str2double(get(handles.peak_detection_threshold,'String'));
settings.min_num_of_peaks=str2double(get(handles.min_num_of_peaks,'String'));
settings.good_fit_cond=str2double(get(handles.good_fit_cond,'String'));
settings.mapping_range_extra=str2double(get(handles.mapping_range_extra,'String'));
settings.autosave_time=str2double(get(handles.autosave_time,'String'));
settings.hmc_threshold=str2double(get(handles.hmc_threshold,'String'));
settings.hmc_num_peaks=str2double(get(handles.hmc_num_peaks,'String'));
settings.min_mol_length=str2double(get(handles.min_mol_length,'String'));
settings.interest_zones=get(handles.interest_zones,'String');
settings.interest_zones_check=logical(get(handles.interest_zones_check,'Value'));
settings.selected_IDs_check=logical(get(handles.selected_IDs_check,'Value'));
settings.selected_IDs_list=get(handles.selected_IDs_list,'String');
settings.align_ch_num=logical(get(handles.align_ch_num,'Value'));
settings.include_non_mapped=logical(get(handles.include_non_mapped,'Value'));
save(get(handles.settings_path, 'String'),'settings');

function background_samples_per_pixel_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function background_samples_per_pixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function background_sub_max_images_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function background_sub_max_images_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function peak_detection_threshold_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function peak_detection_threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_num_of_peaks_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function min_num_of_peaks_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function good_fit_cond_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function good_fit_cond_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function status_bar_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function status_bar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function no_peaks_stat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function no_peaks_stat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function failed_fit_stat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function failed_fit_stat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mapping_OK_stat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function mapping_OK_stat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function range_OK_stat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function range_OK_stat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mapping_range_extra_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function mapping_range_extra_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function autosave_time_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function autosave_time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hmc_marked_stat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function hmc_marked_stat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hmc_threshold_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function hmc_threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in interest_zone_add.
function interest_zone_add_Callback(hObject, eventdata, handles)
interest_zones=get(handles.interest_zones,'String');
zone_to_add=[get(handles.interest_chr,'String'),',',get(handles.interest_start,'String'),',',get(handles.interest_end,'String')];
interest_zones{numel(interest_zones)+1}=zone_to_add;
set(handles.interest_zones,'String',interest_zones);


% --- Executes on selection change in interest_zones.
function interest_zones_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function interest_zones_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in interest_zone_delete.
function interest_zone_delete_Callback(hObject, eventdata, handles)
interest_zones=get(handles.interest_zones,'String');
if length(interest_zones)>1
    interest_zones(get(handles.interest_zones,'value'))=[];
    set(handles.interest_zones,'String',interest_zones,'value',length(interest_zones));
end

function hmc_num_peaks_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function hmc_num_peaks_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_mol_length_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function min_mol_length_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in selected_IDs_check.
function selected_IDs_check_Callback(hObject, eventdata, handles)

function selected_IDs_list_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function selected_IDs_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in writeFOVs.
function writeFOVs_Callback(hObject, eventdata, handles)


% --- Executes on button press in align_ch_num.
function align_ch_num_Callback(hObject, eventdata, handles)
% hObject    handle to align_ch_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of align_ch_num


% --- Executes on button press in include_non_mapped.
function include_non_mapped_Callback(hObject, eventdata, handles)
% hObject    handle to include_non_mapped (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_non_mapped
