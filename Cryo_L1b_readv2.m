function [HDR, CS]=Cryo_L1b_readv2(full_filename)

% Scope: Matlab Function to ingest CryoSat L1b data products to matlab workspace
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original code found: https://earth.esa.int/web/guest/software-tools/
% -/article/cryosat-matlab-routines
% ***** The following Matlab routines read the L1b and L2 CryoSat products.
% ***** The routines have been kindly provided by S. Dinardo, working in 
% ***** the Altimetry Team in ESA/Esrin.
% 
% Update: 2nd October, 2013
% Author: N.Galin (all I did is modify the code to include SARIN FBR
% support).
% Comment: Updated code to include support for SARIN FBR mode data support.
% Supported Modes: FBR (SAR,SARIN), LRM, SAR, SARin, FDM
% Version: 2.0
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supported Modes: FBR, LRM, SAR, SARin, FDM
% 
% Input Argument: <full_filename> is the full pathname (path + file) where the CryoSat .DBL
% file is stored in your local drive
%
% Output Arguments: the structures <HDR>, containing the header of the read file, and the structure <CS>, 
% containing the read data fields
%
% Author: Salvatore Dinardo, Alexander Horvath
% Date: 20/07/2011
% Version: 1.0
% Compliance to CryoSat L1b Product Specification: 4.6 
% Debugging: for any issues, please write to salvatore.dinardo@esa.int

[pathname filename, ext]=fileparts(full_filename);
fid=fopen(fullfile(pathname,[filename ext ] ),'r','b');
s = dir(fullfile(pathname, [filename ext]  ));

%%%%%%%%%%%%%%%%%%%%%%%%  DATA HEADER READING %%%%%%%%%%%%%%%%%%%

MPH_size=1247;

i=0;
k=1;
while 1
    
    i=i+1;
    
    if ftell(fid)>MPH_size, break,   end
    tline = fgetl(fid);
    
    field=tline(1:strfind(tline,'=')-1);
    
    I=strfind(tline,'=')+1;
    
    if  i>2 && isfield(HDR,field)
        
        
        if strcmp(field,'DS_NAME')
            
            k=k+1;
            
        end
        
        field=[field num2str(k)];
        
    end
    
    if strcmp(tline(I),'"')
        
        value=tline(I+1:end-1);
        eval([ 'HDR.' field '='' ' value ''';']);
        
    else
        
        J=strfind(tline,'<')-1;
        
        if isempty(J) 
            J=length(tline);
        end
        
        if  not(isempty(tline(I:J)))&& not(isnan(str2double(tline(I:J))))
            
            value=str2double(tline(I:J));
            eval([ 'HDR.' field '= ' num2str(value, '%10.5f') ';']);
            
        elseif not(isempty(tline(I:J)))&& (isnan(str2double(tline(I:J))))
            
            value=(tline(I:J));
            eval([ 'HDR.' field '= ''' value ''';']);
        end
        
    end
    
end

i=0;
k=1;

while 1
    
    i=i+1;
    
    if ftell(fid)>=MPH_size+HDR.SPH_SIZE 
        break
    end
    
    tline = fgetl(fid);
    field=tline(1:strfind(tline,'=')-1);
    
    I=strfind(tline,'=')+1;
    
    if  i>2 && isfield(HDR,field)
        
        
        if strcmp(field,'DS_NAME')
            
            k=k+1;
            
        end
        
        field=[field num2str(k)];
        
    end
    
    if strcmp(tline(I),'"')
        
        value=tline(I+1:end-1);
        eval([ 'HDR.' field '='' ' value ''';']);
        
    else
        
        J=strfind(tline,'<')-1;
        if isempty(J)
            J=length(tline);
        end
        
        if  not(isempty(tline(I:J)))&& not(isnan(str2double(tline(I:J))))
            
            value=str2double(tline(I:J));
            eval([ 'HDR.' field '= ' num2str(value, '%10.5f') ';']);
            
        elseif not(isempty(tline(I:J)))&& (isnan(str2double(tline(I:J))))
            
            try
                value=(tline(I:J));
                eval([ 'HDR.' field '= ''' value ''';']);
                
            catch
                keyboard
            end
        end
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% OPERATIVE MODE IDENTIFICATION  %%%%%%%%%%%%%%%%%%

CS.GEO.OPERATION_MODE=strtrim(HDR.DS_NAME);
CS.GEO.Start_Time=(datenum(HDR.START_RECORD_TAI_TIME)-datenum('01-Jan-2000 00:00:00')).*24.*60.*60;

%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INFORMATION %%%%%%%%%%%%%%%%%%%%%%

N_block=20;
SAR_pulses_burst=64;
N_samples=128;
N_SIN_samples = 512;
timestamp=4+4+4;
time_group=timestamp+4+2.*2+4.*6+3.*3.*4+4;  % Time and orbit group (84 bytes)
measure_group=8+4.*18+4; % Measurement group (84 bytes)
geo_corr=4.*16;          % Geocorrection group (64 bytes)
offset=geo_corr+(time_group+measure_group)*N_block;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'} 
        
        av_wfm_LRM=(4*6+8+2*N_samples+4*2+2*2); % Average Waveform group (300 bytes )
        wfm_LRM=2*N_samples+4*2+2*2; % Full Waveform group (268 bytes )
        record_size=geo_corr+av_wfm_LRM+(time_group+measure_group+wfm_LRM)*N_block;
        n_points=N_samples;
        n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_FBR_SAR'
        
        wfm_SAR=2*N_samples*SAR_pulses_burst+2+2; % Waveform group (64 bytes )
        record_size=geo_corr+(time_group+measure_group+wfm_SAR)*N_block;
        n_points=2*N_samples*SAR_pulses_burst;
        n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_L1B_SAR'
        
        wfm_SAR=2*N_samples+4*2+2*2+50*2;
        av_wfm_SAR=(timestamp+4*3+8+2*N_samples+4*2+2*2);
        record_size=geo_corr+av_wfm_SAR+(time_group+measure_group+wfm_SAR)*N_block;
        n_points=N_samples;
        n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
   
    case 'SIR_FBR_SARIN'
        
        wfm_SIN=2*N_SIN_samples*SAR_pulses_burst+2*N_SIN_samples*SAR_pulses_burst+2+2; % Waveform group (64 bytes ) %TODO
        record_size=geo_corr+(time_group+measure_group+wfm_SIN)*N_block; 
        n_points=2*N_SIN_samples*SAR_pulses_burst+2*N_SIN_samples*SAR_pulses_burst;
        n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_L1B_SARIN'
    
        wfm_SIN=512*2+4+4+2+2+50*2+N_samples.*4.*2+N_samples.*4.*4;
        av_wfm_SIN=12+4+4+4+8+N_samples.*4.*2+4+4+2+2;
        record_size=geo_corr+av_wfm_SIN+(time_group+measure_group+wfm_SIN)*N_block;
        n_points=N_samples.*4; % 512 bins insted of 128 bins
        n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
end


%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%%%

CS.GEO.TAI.days=zeros(N_block,n_recs);
CS.GEO.TAI.secs=zeros(N_block,n_recs);
CS.GEO.TAI.microsecs=zeros(N_block,n_recs);
CS.GEO.USO=zeros(N_block,n_recs);
CS.GEO.MODE_ID=char(zeros(N_block,16,n_recs));
CS.GEO.SRC_CNT=zeros(N_block,n_recs);
CS.GEO.INS_CFG=char(zeros(N_block,32,n_recs));
CS.GEO.BURST_CNT=zeros(N_block,n_recs);
CS.GEO.LAT=zeros(N_block,n_recs);
CS.GEO.LON=zeros(N_block,n_recs);
CS.GEO.H=zeros(N_block,n_recs);
CS.GEO.H_rate=zeros(N_block,n_recs);
CS.GEO.V.Vx=zeros(N_block,n_recs);
CS.GEO.V.Vy=zeros(N_block,n_recs);
CS.GEO.V.Vz=zeros(N_block,n_recs);
CS.GEO.Beam.X=zeros(N_block,n_recs);
CS.GEO.Beam.Y=zeros(N_block,n_recs);
CS.GEO.Beam.Z=zeros(N_block,n_recs);
CS.GEO.BaseLine.X=zeros(N_block,n_recs);
CS.GEO.BaseLine.Y=zeros(N_block,n_recs);
CS.GEO.BaseLine.Z=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG=char(zeros(N_block,32,n_recs));



CS.GEO.MODE_ID_Tab.Ins_Oper_Mode=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Sarin_Degr_Case=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Reserved_1=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Cal4_Flag=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Plat_Att_Ctrl=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Reserved_2=zeros(N_block,n_recs);


CS.GEO.INS_CFG_Tab.Rx_Chain_Use=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.SIRAL_ID=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_1=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Band_FLAG=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_2=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Tracking_Mode=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Ext_Cal=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_3=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Loop_Status=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Echo_Loss=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Real_Time_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Echo_Satur_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Rx_Band_Atten=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Cycle_Report_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Star_Trk_1=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Star_Trk_2=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Star_Trk_3=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_4=zeros(N_block,n_recs);


CS.GEO.MCD_FLAG_Tab.Block_Degraded=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Blank_Block=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Datation_Degraded=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Orbit_Propag_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Orbit_File_Change=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Orbit_Discontinuity=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Echo_Saturation=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Other_Echo_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Rx1_Err_SARin=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Rx2_Err_SARin=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Wind_Delay_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.AGC_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.CAL1_Corr_Miss=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.CAL1_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.DORIS_USO_Corr=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Complex_CAL1_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.TRK_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.RX1_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.RX2_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.NPM_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Reserved_1=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Reserved_2=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Reserved_3=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Reserved_4=zeros(N_block,n_recs);

CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Corr=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.CAL2_Corr_Miss=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.CAL2_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Power_Scaling_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Att_Corr_Miss=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Reserved_5=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Reserved_6=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Mode=zeros(N_block,n_recs);

CS.MEA.win_delay=zeros(N_block,n_recs);
CS.MEA.Ho=zeros(N_block,n_recs);
CS.MEA.Trk_H_rate=zeros(N_block,n_recs);
CS.MEA.LAI=zeros(N_block,n_recs);
CS.MEA.FAI=zeros(N_block,n_recs);
CS.MEA.AGC_1=zeros(N_block,n_recs);
CS.MEA.AGC_2=zeros(N_block,n_recs);
CS.MEA.Gain_Rx1=zeros(N_block,n_recs);
CS.MEA.Gain_Rx2=zeros(N_block,n_recs);
CS.MEA.Tx_Pwr=zeros(N_block,n_recs);
CS.MEA.dpl_range_corr=zeros(N_block,n_recs);
CS.MEA.ins_range_corr_rx_tx=zeros(N_block,n_recs);
CS.MEA.ins_range_corr_rx=zeros(N_block,n_recs);
CS.MEA.ins_gain_corr_rx_tx=zeros(N_block,n_recs);
CS.MEA.ins_gain_corr_rx=zeros(N_block,n_recs);
CS.MEA.int_phase_corr=zeros(N_block,n_recs);
CS.MEA.ext_phase_corr=zeros(N_block,n_recs);
CS.MEA.noise_power=zeros(N_block,n_recs);
CS.MEA.phase_slope_corr=zeros(N_block,n_recs);

CS.COR.dry_trop=zeros(1,n_recs);
CS.COR.wet_trop=zeros(1,n_recs);
CS.COR.inv_bar=zeros(1,n_recs);
CS.COR.dac=zeros(1,n_recs);
CS.COR.gim_ion=zeros(1,n_recs);
CS.COR.model_ion=zeros(1,n_recs);
CS.COR.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.ocean_loading_tide=zeros(1,n_recs);
CS.COR.solidearth_tide=zeros(1,n_recs);
CS.COR.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.surf_type=zeros(1,n_recs);
CS.COR.corr_status=char(zeros(32,n_recs));
CS.COR.corr_error=char(zeros(32,n_recs));


CS.COR.corr_status_tab.dry_trop=zeros(1,n_recs);
CS.COR.corr_status_tab.wet_trop=zeros(1,n_recs);
CS.COR.corr_status_tab.inv_bar=zeros(1,n_recs);
CS.COR.corr_status_tab.dac=zeros(1,n_recs);
CS.COR.corr_status_tab.gim_iono=zeros(1,n_recs);
CS.COR.corr_status_tab.model_iono=zeros(1,n_recs);
CS.COR.corr_status_tab.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.ocean_loading_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.solidearth_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.surface_type=zeros(1,n_recs);
CS.COR.corr_status_tab.reserved=zeros(1,n_recs);

CS.COR.corr_error_tab.dry_trop=zeros(1,n_recs);
CS.COR.corr_error_tab.wet_trop=zeros(1,n_recs);
CS.COR.corr_error_tab.inv_bar=zeros(1,n_recs);
CS.COR.corr_error_tab.dac=zeros(1,n_recs);
CS.COR.corr_error_tab.gim_iono=zeros(1,n_recs);
CS.COR.corr_error_tab.model_iono=zeros(1,n_recs);
CS.COR.corr_error_tab.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.ocean_loading_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.solidearth_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.surface_type=zeros(1,n_recs);
CS.COR.corr_error_tab.reserved=zeros(1,n_recs);


switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'} 
        
        
        CS.AVG.TAI.day=zeros(1,n_recs);
        CS.AVG.TAI.secs=zeros(1,n_recs);
        CS.AVG.TAI.microsecs=zeros(1,n_recs);
        CS.AVG.lat=zeros(1,n_recs);
        CS.AVG.lon=zeros(1,n_recs);
        CS.AVG.H=zeros(1,n_recs);
        CS.AVG.win_delay=zeros(1,n_recs);
        CS.AVG.data=zeros(n_points,n_recs);
        CS.AVG.echo_scaling=zeros(1,n_recs);
        CS.AVG.echo_scale_power=zeros(1,n_recs);
        CS.AVG.N_averaged_echoes=zeros(1,n_recs);
        CS.AVG.OneHz_Echo_Err=zeros(1,n_recs);
        CS.LRM.data=zeros(n_points.*N_block,n_recs);
        CS.LRM.echo_scaling=zeros(N_block,n_recs);
        CS.LRM.echo_scale_power=zeros(N_block,n_recs);
        CS.LRM.N_averaged_echoes=zeros(N_block,n_recs);
        CS.LRM.FLAG=char(zeros(16,N_block,n_recs));

    case {'SIR_L1B_SAR'}
        
        
        CS.AVG.TAI.days=zeros(1,n_recs);
        CS.AVG.TAI.secs=zeros(1,n_recs);
        CS.AVG.TAI.microsecs=zeros(1,n_recs);
        CS.AVG.lat=zeros(1,n_recs);
        CS.AVG.lon=zeros(1,n_recs);
        CS.AVG.H=zeros(1,n_recs);
        CS.AVG.win_delay=zeros(1,n_recs);
        CS.AVG.data=zeros(n_points,n_recs);
        CS.AVG.echo_scaling=zeros(1,n_recs);
        CS.AVG.echo_scale_power=zeros(1,n_recs);
        CS.AVG.N_averaged_echoes=zeros(1,n_recs);
        CS.AVG.OneHz_Echo_Err=zeros(1,n_recs);
        CS.AVG.Mispointing_Err=zeros(1,n_recs);

        CS.SAR.data=zeros(n_points.*N_block,n_recs);
        CS.SAR.echo_scaling=zeros(N_block,n_recs);
        CS.SAR.echo_scale_power=zeros(N_block,n_recs);
        CS.SAR.N_averaged_echoes=zeros(N_block,n_recs);
        CS.SAR.FLAG=char(zeros(16,N_block,n_recs));
        CS.SAR.FLAG_tab.Approximate_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Exact_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Doppler_Weighting_Computed=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Doppler_Weighting_Applied_Before_Stack=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Multilook_Incomplete=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Beam_Angle_Steering_Err=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.AntiAliased_Power_Echo=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Auto_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Reserved=zeros(N_block,n_recs);
        CS.SAR.beam_param=zeros(50.*N_block,n_recs);
        
        
    case {'SIR_FBR_SAR'}
        
        CS.FBR.data=zeros(n_points.*N_block,n_recs,'int8');
        CS.FBR.N_pulses=zeros(N_block,n_recs);
        CS.FBR.FLAG=char(zeros(16,N_block,n_recs));
        
    case {'SIR_FBR_SARIN'}
        
        CS.FBR.data=zeros(n_points.*N_block,n_recs,'int8');
        CS.FBR.Rx1=zeros(N_SIN_samples,N_block,SAR_pulses_burst,n_recs,'int8');
        CS.FBR.Rx2=zeros(N_SIN_samples,N_block,SAR_pulses_burst,n_recs,'int8');
        CS.FBR.N_pulses=zeros(N_block,n_recs);
        CS.FBR.FLAG=char(zeros(16,N_block,n_recs));
        
    case 'SIR_L1B_SARIN'
        
        CS.AVG.TAI.days=zeros(1,n_recs);
        CS.AVG.TAI.secs=zeros(1,n_recs);
        CS.AVG.TAI.microsecs=zeros(1,n_recs);
        CS.AVG.lat=zeros(1,n_recs);
        CS.AVG.lon=zeros(1,n_recs);
        CS.AVG.H=zeros(1,n_recs);
        CS.AVG.win_delay=zeros(1,n_recs);
        CS.AVG.data=zeros(n_points,n_recs); % n_points in this case 512 instead of 128
        CS.AVG.echo_scaling=zeros(1,n_recs);
        CS.AVG.echo_scale_power=zeros(1,n_recs);
        CS.AVG.N_averaged_echoes=zeros(1,n_recs);
        CS.AVG.OneHz_Echo_Err=zeros(1,n_recs);
        CS.AVG.Mispointing_Err=zeros(1,n_recs);

        CS.SIN.data=zeros(n_points.*N_block,n_recs); % n_points in this case 512 instead of 128
        CS.SIN.echo_scaling=zeros(N_block,n_recs);
        CS.SIN.echo_scale_power=zeros(N_block,n_recs);
        CS.SIN.N_averaged_echoes=zeros(N_block,n_recs);
        CS.SIN.FLAG=char(zeros(16,N_block,n_recs));
        CS.SIN.FLAG_tab.Approximate_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Exact_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Doppler_Weighting_Computed=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Doppler_Weighting_Applied_Before_Stack=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Multilook_Incomplete=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Beam_Angle_Steering_Err=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.AntiAliased_Power_Echo=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Auto_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Reserved=zeros(N_block,n_recs);
        CS.SIN.beam_param=zeros(50.*N_block,n_recs);
        CS.SIN.coherence=zeros(n_points.*N_block,n_recs);
        CS.SIN.phase_difference=zeros(n_points.*N_block,n_recs);
        
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:n_recs
    
    
    %%%%%%%%% Time and Orbit Group Reading %%%%%%%%%%%%%%%%%%%%%%%%
    
    fseek(fid,MPH_size+HDR.SPH_SIZE+(i-1).*record_size,'bof');
    CS.GEO.TAI.days(:,i)=fread(fid,N_block,'int32',time_group-4);
    fseek(fid,MPH_size+HDR.SPH_SIZE+4+(i-1).*record_size,'bof');
    CS.GEO.TAI.secs(:,i)=fread(fid,N_block,'uint32',time_group-4);
    fseek(fid,MPH_size+HDR.SPH_SIZE+4+4+(i-1).*record_size,'bof');
    CS.GEO.TAI.microsecs(:,i)=fread(fid,N_block,'uint32',time_group-4);
    fseek(fid,MPH_size+HDR.SPH_SIZE+timestamp+(i-1).*record_size,'bof');
    CS.GEO.USO(:,i)=fread(fid,N_block,'int32',time_group-4).*1e-9;
    fseek(fid,MPH_size+HDR.SPH_SIZE+timestamp+ 4*1+(i-1).*record_size,'bof');
    CS.GEO.MODE_ID(:,:,i)=dec2bin(fread(fid,N_block ,'uint16',time_group-2),16);
    fseek(fid,MPH_size+HDR.SPH_SIZE+timestamp+ 4*1+2+(i-1).*record_size,'bof');
    CS.GEO.SRC_CNT(:,i)=fread(fid,N_block ,'uint16',time_group-2)+1;
    fseek(fid,MPH_size+HDR.SPH_SIZE+timestamp+ 4*1+2+2+(i-1).*record_size,'bof');
    CS.GEO.INS_CFG(:,:,i)=dec2bin(fread(fid,N_block ,'uint32',time_group-4),32);
    fseek(fid,MPH_size+HDR.SPH_SIZE+timestamp+ 4*2+2+2+(i-1).*record_size,'bof');
    CS.GEO.BURST_CNT(:,i)=fread(fid,N_block ,'uint32',time_group-4);
    fseek(fid, MPH_size+HDR.SPH_SIZE+timestamp+4*3+2*2 +(i-1).*record_size,'bof');
    CS.GEO.LAT(:,i)=fread(fid,N_block,'int32',time_group-4)./1e7;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*4+2*2 +(i-1).*record_size,'bof');
    CS.GEO.LON(:,i)=fread(fid,N_block,'int32',time_group-4)./1e7;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*5+2*2 +(i-1).*record_size,'bof');
    CS.GEO.H(:,i)=fread(fid,N_block,'int32',time_group-4)./1e3;
    fseek(fid, MPH_size+HDR.SPH_SIZE+timestamp+4*6+2*2 +(i-1).*record_size,'bof');
    CS.GEO.H_rate(:,i)=fread(fid,N_block,'int32',time_group-4)./1e3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*7+2*2 +(i-1).*record_size,'bof');
    CS.GEO.V.Vx(:,i)=fread(fid,N_block,'int32',time_group-4)./1e3;
    fseek(fid, MPH_size+HDR.SPH_SIZE+timestamp+4*8+2*2 +(i-1).*record_size,'bof');
    CS.GEO.V.Vy(:,i)=fread(fid,N_block,'int32',time_group-4)./1e3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*9+2*2 +(i-1).*record_size,'bof');
    CS.GEO.V.Vz(:,i)=fread(fid,N_block,'int32',time_group-4)./1e3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*10+2*2 +(i-1).*record_size,'bof');
    CS.GEO.Beam.X(:,i)=fread(fid,N_block,'int32',time_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*11+2*2 +(i-1).*record_size,'bof');
    CS.GEO.Beam.Y(:,i)=fread(fid,N_block,'int32',time_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*12+2*2 +(i-1).*record_size,'bof');
    CS.GEO.Beam.Z(:,i)=fread(fid,N_block,'int32',time_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*13+2*2 +(i-1).*record_size,'bof');
    CS.GEO.BaseLine.X(:,i)=fread(fid,N_block,'int32',time_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+timestamp+4*14+2*2 +(i-1).*record_size,'bof');
    CS.GEO.BaseLine.Y(:,i)=fread(fid,N_block,'int32',time_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*15+2*2 +(i-1).*record_size,'bof');
    CS.GEO.BaseLine.Z(:,i)=fread(fid,N_block,'int32',time_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+ timestamp+4*16+2*2 +(i-1).*record_size,'bof');
    CS.GEO.MCD_FLAG(:,:,i)=dec2bin(fread(fid,N_block,'uint32',time_group-4),32);

    CS.GEO.MODE_ID_Tab.Ins_Oper_Mode(:,i)=bin2dec(CS.GEO.MODE_ID(:,1:6,i));
    CS.GEO.MODE_ID_Tab.Sarin_Degr_Case(:,i)=bin2dec(CS.GEO.MODE_ID(:,7,i));
    CS.GEO.MODE_ID_Tab.Reserved_1(:,i)=bin2dec(CS.GEO.MODE_ID(:,8,i));
    CS.GEO.MODE_ID_Tab.Cal4_Flag(:,i)=bin2dec(CS.GEO.MODE_ID(:,9,i));
    CS.GEO.MODE_ID_Tab.Plat_Att_Ctrl(:,i)=bin2dec(CS.GEO.MODE_ID(:,10:11,i));
    CS.GEO.MODE_ID_Tab.Reserved_2(:,i)=bin2dec(CS.GEO.MODE_ID(:,12:16,i));

    CS.GEO.INS_CFG_Tab.Rx_Chain_Use(:,i)=bin2dec(CS.GEO.INS_CFG(:,1:2,i));
    CS.GEO.INS_CFG_Tab.SIRAL_ID(:,i)=bin2dec(CS.GEO.INS_CFG(:,3,i));
    CS.GEO.INS_CFG_Tab.Reserved_1(:,i)=bin2dec(CS.GEO.INS_CFG(:,4,i));
    CS.GEO.INS_CFG_Tab.Band_FLAG(:,i)=bin2dec(CS.GEO.INS_CFG(:,5:6,i));
    CS.GEO.INS_CFG_Tab.Reserved_2(:,i)=bin2dec(CS.GEO.INS_CFG(:,7:8,i));
    CS.GEO.INS_CFG_Tab.Tracking_Mode(:,i)=bin2dec(CS.GEO.INS_CFG(:,9:10,i));
    CS.GEO.INS_CFG_Tab.Ext_Cal(:,i)=bin2dec(CS.GEO.INS_CFG(:,11,i));
    CS.GEO.INS_CFG_Tab.Reserved_3(:,i)=bin2dec(CS.GEO.INS_CFG(:,12,i));
    CS.GEO.INS_CFG_Tab.Loop_Status(:,i)=bin2dec(CS.GEO.INS_CFG(:,13,i));
    CS.GEO.INS_CFG_Tab.Echo_Loss(:,i)=bin2dec(CS.GEO.INS_CFG(:,14,i));
    CS.GEO.INS_CFG_Tab.Real_Time_Err(:,i)=bin2dec(CS.GEO.INS_CFG(:,15,i));
    CS.GEO.INS_CFG_Tab.Echo_Satur_Err(:,i)=bin2dec(CS.GEO.INS_CFG(:,16,i));
    CS.GEO.INS_CFG_Tab.Rx_Band_Atten(:,i)=bin2dec(CS.GEO.INS_CFG(:,17,i));
    CS.GEO.INS_CFG_Tab.Cycle_Report_Err(:,i)=bin2dec(CS.GEO.INS_CFG(:,18,i));
    CS.GEO.INS_CFG_Tab.Star_Trk_1(:,i)=bin2dec(CS.GEO.INS_CFG(:,19,i));
    CS.GEO.INS_CFG_Tab.Star_Trk_2(:,i)=bin2dec(CS.GEO.INS_CFG(:,20,i));
    CS.GEO.INS_CFG_Tab.Star_Trk_3(:,i)=bin2dec(CS.GEO.INS_CFG(:,21,i));
    CS.GEO.INS_CFG_Tab.Reserved_4(:,i)=bin2dec(CS.GEO.INS_CFG(:,22:32,i));

    
    CS.GEO.MCD_FLAG_Tab.Block_Degraded(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,1,i));
    CS.GEO.MCD_FLAG_Tab.Blank_Block(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,2,i));
    CS.GEO.MCD_FLAG_Tab.Datation_Degraded(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,3,i));
    CS.GEO.MCD_FLAG_Tab.Orbit_Propag_Err(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,4,i));
    CS.GEO.MCD_FLAG_Tab.Orbit_File_Change(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,5,i));
    CS.GEO.MCD_FLAG_Tab.Orbit_Discontinuity(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,6,i));
    CS.GEO.MCD_FLAG_Tab.Echo_Saturation(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,7,i));
    CS.GEO.MCD_FLAG_Tab.Other_Echo_Err(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,8,i));
    CS.GEO.MCD_FLAG_Tab.Rx1_Err_SARin(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,9,i));
    CS.GEO.MCD_FLAG_Tab.Rx2_Err_SARin(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,10,i));
    CS.GEO.MCD_FLAG_Tab.Wind_Delay_Incon(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,11,i));
    CS.GEO.MCD_FLAG_Tab.AGC_Incon(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,12,i));
    CS.GEO.MCD_FLAG_Tab.CAL1_Corr_Miss(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,13,i));
    CS.GEO.MCD_FLAG_Tab.CAL1_Corr_IPF(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,14,i));
    CS.GEO.MCD_FLAG_Tab.DORIS_USO_Corr(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,15,i));
    CS.GEO.MCD_FLAG_Tab.Complex_CAL1_Corr_IPF(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,16,i));
    CS.GEO.MCD_FLAG_Tab.TRK_ECHO_Err(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,17,i));
    CS.GEO.MCD_FLAG_Tab.RX1_ECHO_Err(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,18,i));
    CS.GEO.MCD_FLAG_Tab.RX2_ECHO_Err(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,19,i));
    CS.GEO.MCD_FLAG_Tab.NPM_Incon(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,20,i));
    CS.GEO.MCD_FLAG_Tab.Reserved_1(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,21,i));
    CS.GEO.MCD_FLAG_Tab.Reserved_2(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,22,i));
    CS.GEO.MCD_FLAG_Tab.Reserved_3(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,23,i));
    CS.GEO.MCD_FLAG_Tab.Reserved_4(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,24,i));
    CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Corr(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,25,i));
    CS.GEO.MCD_FLAG_Tab.CAL2_Corr_Miss(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,26,i));
    CS.GEO.MCD_FLAG_Tab.CAL2_Corr_IPF(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,27,i));
    CS.GEO.MCD_FLAG_Tab.Power_Scaling_Err(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,28,i));
    CS.GEO.MCD_FLAG_Tab.Att_Corr_Miss(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,29,i));
    CS.GEO.MCD_FLAG_Tab.Reserved_5(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,30,i));  
    CS.GEO.MCD_FLAG_Tab.Reserved_6(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,31,i));  
    CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Mode(:,i)=bin2dec(CS.GEO.MCD_FLAG(:,32,i));  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%% Measurements Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+(i-1).*record_size,'bof');
    CS.MEA.win_delay(:,i)=fread(fid,N_block,'int64',measure_group-8).*1e-12;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+(i-1).*record_size,'bof');
    CS.MEA.Ho(:,i)=fread(fid,N_block,'int32',measure_group-4).*48.8e-12;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*1+(i-1).*record_size,'bof');
    CS.MEA.Trk_H_rate(:,i)=fread(fid,N_block,'int32',measure_group-4);
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*2+(i-1).*record_size,'bof');
    CS.MEA.LAI(:,i)=fread(fid,N_block,'int32',measure_group-4).*12.5e-9;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*3+(i-1).*record_size,'bof');
    CS.MEA.FAI(:,i)=fread(fid,N_block,'int32',measure_group-4).*12.5./256.*1e-9;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*4+(i-1).*record_size,'bof');
    CS.MEA.AGC_1(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*5+(i-1).*record_size,'bof');
    CS.MEA.AGC_2(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*6+(i-1).*record_size,'bof');
    CS.MEA.Gain_Rx1(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*7+(i-1).*record_size,'bof');
    CS.MEA.Gain_Rx2(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*8+(i-1).*record_size,'bof');
    CS.MEA.Tx_Pwr(:,i)=fread(fid,N_block,'int32',measure_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*9+(i-1).*record_size,'bof');
    CS.MEA.dpl_range_corr(:,i)=fread(fid,N_block,'int32',measure_group-4)./1000;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*10+(i-1).*record_size,'bof');
    CS.MEA.ins_range_corr_rx_tx(:,i)=fread(fid,N_block,'int32',measure_group-4)./1000;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*11+(i-1).*record_size,'bof');
    CS.MEA.ins_range_corr_rx(:,i)=fread(fid,N_block,'int32',measure_group-4)./1000;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*12+(i-1).*record_size,'bof');
    CS.MEA.ins_gain_corr_rx_tx(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*13+(i-1).*record_size,'bof');
    CS.MEA.ins_gain_corr_rx(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*14+(i-1).*record_size,'bof');
    CS.MEA.int_phase_corr(:,i)=fread(fid,N_block,'int32',measure_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*15+(i-1).*record_size,'bof');
    CS.MEA.ext_phase_corr(:,i)=fread(fid,N_block,'int32',measure_group-4)./1e6;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*16+(i-1).*record_size,'bof');
    CS.MEA.noise_power(:,i)=fread(fid,N_block,'int32',measure_group-4)./100;
    fseek(fid,MPH_size+HDR.SPH_SIZE+time_group.*N_block+8+4*17+(i-1).*record_size,'bof');
    CS.MEA.phase_slope_corr(:,i)=fread(fid,N_block,'int32',measure_group-4)./1e6;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%% Corrections Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4.*0+(i-1).*record_size,'bof');
    CS.COR.dry_trop(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*1+(i-1).*record_size,'bof');
    CS.COR.wet_trop(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*2+(i-1).*record_size,'bof');
    CS.COR.inv_bar(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*3+(i-1).*record_size,'bof');
    CS.COR.dac(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*4+(i-1).*record_size,'bof');
    CS.COR.gim_ion(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*5+(i-1).*record_size,'bof');
    CS.COR.model_ion(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*6+(i-1).*record_size,'bof');
    CS.COR.ocean_equilibrium_tide(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*7+(i-1).*record_size,'bof');
    CS.COR.ocean_longperiod_tide(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*8+(i-1).*record_size,'bof');
    CS.COR.ocean_loading_tide(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*9+(i-1).*record_size,'bof');
    CS.COR.solidearth_tide(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*10+(i-1).*record_size,'bof');
    CS.COR.geocentric_polar_tide(i)=fread(fid,1,'int32').*1e-3;
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*11+(i-1).*record_size,'bof');
    CS.COR.surf_type(i)=fread(fid,1,'uint32');
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*13+(i-1).*record_size,'bof');
    CS.COR.corr_status(:,i)=dec2bin(fread(fid,1,'uint32'),32).';
    fseek(fid,MPH_size+HDR.SPH_SIZE+(time_group+measure_group)*N_block+4*14+(i-1).*record_size,'bof');
    CS.COR.corr_error(:,i)=dec2bin(fread(fid,1,'uint32'),32).';
    
    
   CS.COR.corr_status_tab.dry_trop(i)=bin2dec(CS.COR.corr_status(1,i));
   CS.COR.corr_status_tab.wet_trop(i)=bin2dec(CS.COR.corr_status(2,i));
   CS.COR.corr_status_tab.inv_bar(i)=bin2dec(CS.COR.corr_status(3,i));
   CS.COR.corr_status_tab.dac(i)=bin2dec(CS.COR.corr_status(4,i));
   CS.COR.corr_status_tab.gim_iono(i)=bin2dec(CS.COR.corr_status(5,i));
   CS.COR.corr_status_tab.model_iono(i)=bin2dec(CS.COR.corr_status(6,i));
   CS.COR.corr_status_tab.ocean_equilibrium_tide(i)=bin2dec(CS.COR.corr_status(7,i));
   CS.COR.corr_status_tab.ocean_longperiod_tide(i)=bin2dec(CS.COR.corr_status(8,i));
   CS.COR.corr_status_tab.ocean_loading_tide(i)=bin2dec(CS.COR.corr_status(9,i));
   CS.COR.corr_status_tab.solidearth_tide(i)=bin2dec(CS.COR.corr_status(10,i));
   CS.COR.corr_status_tab.geocentric_polar_tide(:,i)=bin2dec(CS.COR.corr_status(11,i));
   CS.COR.corr_status_tab.surface_type(i)=bin2dec(CS.COR.corr_status(12,i));
   CS.COR.corr_status_tab.reserved(i)=bin2dec(CS.COR.corr_status(13:32,i).');
    
    
   CS.COR.corr_error_tab.dry_trop(i)=bin2dec(CS.COR.corr_error(1,i));
   CS.COR.corr_error_tab.wet_trop(i)=bin2dec(CS.COR.corr_error(2,i));
   CS.COR.corr_error_tab.inv_bar(i)=bin2dec(CS.COR.corr_error(3,i));
   CS.COR.corr_error_tab.dac(i)=bin2dec(CS.COR.corr_error(4,i));
   CS.COR.corr_error_tab.gim_iono(i)=bin2dec(CS.COR.corr_error(5,i));
   CS.COR.corr_error_tab.model_iono(i)=bin2dec(CS.COR.corr_error(6,i));
   CS.COR.corr_error_tab.ocean_equilibrium_tide(i)=bin2dec(CS.COR.corr_error(7,i));
   CS.COR.corr_error_tab.ocean_longperiod_tide(i)=bin2dec(CS.COR.corr_error(8,i));
   CS.COR.corr_error_tab.ocean_loading_tide(i)=bin2dec(CS.COR.corr_error(9,i));
   CS.COR.corr_error_tab.solidearth_tide(i)=bin2dec(CS.COR.corr_error(10,i));
   CS.COR.corr_error_tab.geocentric_polar_tide(:,i)=bin2dec(CS.COR.corr_error(11,i));
   CS.COR.corr_error_tab.surface_type(i)=bin2dec(CS.COR.corr_error(12,i));
   CS.COR.corr_error_tab.reserved(i)=bin2dec(CS.COR.corr_error(13:32,i).');
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    switch CS.GEO.OPERATION_MODE
        
        case {'SIR_L1B_LRM','SIR_L1B_FDM'} 
            
            %%%%%%%%% LRM Average Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*0+(i-1).*record_size,'bof');
            CS.AVG.TAI.day(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*1+(i-1).*record_size,'bof');
            CS.AVG.TAI.secs(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*2+(i-1).*record_size,'bof');
            CS.AVG.TAI.microsecs(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*3+(i-1).*record_size,'bof');
            CS.AVG.lat(:,i)=fread(fid,1,'int32').*1e-7;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*4+(i-1).*record_size,'bof');
            CS.AVG.lon(:,i)=fread(fid,1,'int32').*1e-7;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*5+(i-1).*record_size,'bof');
            CS.AVG.H(:,i)=fread(fid,1,'int32')*1e-3;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+(i-1).*record_size,'bof');
            CS.AVG.win_delay(:,i)=fread(fid,1,'uint64')*1e-12;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+(i-1).*record_size,'bof');
            CS.AVG.data(:,i)=fread(fid,n_points,'uint16');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+2.*n_points+(i-1).*record_size,'bof');
            CS.AVG.echo_scaling(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+2.*n_points+4.*1+(i-1).*record_size,'bof');
            CS.AVG.echo_scale_power(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset++4*6+8+2.*n_points+4.*2+(i-1).*record_size,'bof');
            CS.AVG.N_averaged_echoes(:,i)=fread(fid,1,'int16');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4*6+8+2.*n_points+4.*2+2.*1+(i-1).*record_size,'bof');
            dummy=dec2bin(fread(fid,1,'uint16').',16);
            CS.AVG.OneHz_Echo_Err(:,i)=bin2dec(dummy(1).');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%% LRM Full Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_LRM+(i-1).*record_size,'bof');
            CS.LRM.data(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*uint16'],4.*2+2.*2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_LRM+2.*n_points+(i-1).*record_size,'bof');
            CS.LRM.echo_scaling(:,i)=fread(fid,N_block,'int32',wfm_LRM-4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_LRM+2.*n_points+4.*1+(i-1).*record_size,'bof');
            CS.LRM.echo_scale_power(:,i)=fread(fid,N_block,'int32',wfm_LRM-4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_LRM+2.*n_points+4*2+(i-1).*record_size,'bof');
            CS.LRM.N_averaged_echoes(:,i)=fread(fid,N_block,'int16',wfm_LRM-2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_LRM+2.*n_points+4.*2+2.*1+(i-1).*record_size,'bof');
            CS.LRM.FLAG(:,:,i)=dec2bin(fread(fid,N_block,'uint16',wfm_LRM-2).',16).';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case {'SIR_L1B_SAR'}
            
            
            %%%%%%%%% SAR Average Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*0+(i-1).*record_size,'bof');
            CS.AVG.TAI.days(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*1+(i-1).*record_size,'bof');
            CS.AVG.TAI.secs(:,i)=fread(fid,1,'uint32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*2+(i-1).*record_size,'bof');
            CS.AVG.TAI.microsecs(:,i)=fread(fid,1,'uint32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*3+(i-1).*record_size,'bof');
            CS.AVG.lat(:,i)=fread(fid,1,'int32').*1e-7;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*4+(i-1).*record_size,'bof');
            CS.AVG.lon(:,i)=fread(fid,1,'int32').*1e-7;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*5+(i-1).*record_size,'bof');
            CS.AVG.H(:,i)=fread(fid,1,'int32')*1e-3;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+(i-1).*record_size,'bof');
            CS.AVG.win_delay(:,i)=fread(fid,1,'uint64')*1e-12;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+(i-1).*record_size,'bof');
            CS.AVG.data(:,i)=fread(fid,n_points,'uint16');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+2.*n_points+(i-1).*record_size,'bof');
            CS.AVG.echo_scaling(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+2.*n_points+4.*1+(i-1).*record_size,'bof');
            CS.AVG.echo_scale_power(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4*6+8+2.*n_points+4.*2+(i-1).*record_size,'bof');
            CS.AVG.N_averaged_echoes(:,i)=fread(fid,1,'uint16');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4*6+8+2.*n_points+4.*2+2.*1+(i-1).*record_size,'bof');
            dummy=dec2bin(fread(fid,1,'uint16').',16);
            CS.AVG.OneHz_Echo_Err(:,i)=bin2dec(dummy(1).');
            CS.AVG.Mispointing_Err(:,i)=bin2dec(dummy(16).');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%% SAR Full Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SAR+(i-1).*record_size,'bof');
            CS.SAR.data(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*uint16'],4.*2+2.*2+50*2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SAR+2.*n_points+(i-1).*record_size,'bof');
            CS.SAR.echo_scaling(:,i)=fread(fid,N_block,'int32',wfm_SAR-4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SAR+2.*n_points+4.*1+(i-1).*record_size,'bof');
            CS.SAR.echo_scale_power(:,i)=fread(fid,N_block,'int32',wfm_SAR-4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SAR+2.*n_points+4*2+(i-1).*record_size,'bof');
            CS.SAR.N_averaged_echoes(:,i)=fread(fid,N_block,'int16',wfm_SAR-2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SAR+2.*n_points+4.*2+2.*1+(i-1).*record_size,'bof');
            CS.SAR.FLAG(:,:,i)=dec2bin(fread(fid,N_block,'uint16',wfm_SAR-2).',16).';
            CS.SAR.FLAG_tab.Approximate_Beam_Steering(:,i)=bin2dec(CS.SAR.FLAG(1,:,i).');
            CS.SAR.FLAG_tab.Exact_Beam_Steering(:,i)=bin2dec(CS.SAR.FLAG(2,:,i).');
            CS.SAR.FLAG_tab.Doppler_Weighting_Computed(:,i)=bin2dec(CS.SAR.FLAG(3,:,i).');
            CS.SAR.FLAG_tab.Doppler_Weighting_Applied_Before_Stack(:,i)=bin2dec(CS.SAR.FLAG(4,:,i).');
            CS.SAR.FLAG_tab.Multilook_Incomplete(:,i)=bin2dec(CS.SAR.FLAG(5,:,i).');
            CS.SAR.FLAG_tab.Beam_Angle_Steering_Err(:,i)=bin2dec(CS.SAR.FLAG(6,:,i).');
            CS.SAR.FLAG_tab.AntiAliased_Power_Echo(:,i)=bin2dec(CS.SAR.FLAG(7,:,i).');
            CS.SAR.FLAG_tab.Auto_Beam_Steering(:,i)=bin2dec(CS.SAR.FLAG(8,:,i).');
            CS.SAR.FLAG_tab.Reserved(:,i)=bin2dec(CS.SAR.FLAG(9:16,:,i).');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SAR+2.*n_points+4*2+2.*2+(i-1).*record_size,'bof');
            CS.SAR.beam_param(:,i)=fread(fid,50.*N_block,[num2str(50) '*uint16'],wfm_SAR-2.*50);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        case {'SIR_FBR_SAR'}
            
            
            %%%%%%%%% FBR Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+(i-1).*record_size,'bof');
            CS.FBR.data(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*int8=>int8'],4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+n_points+(i-1).*record_size,'bof');
            CS.FBR.N_pulses(:,i)=fread(fid,N_block,'uint16',wfm_SAR-2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+n_points+2+(i-1).*record_size,'bof');
            CS.FBR.FLAG(:,:,i)=dec2bin(fread(fid,N_block,'uint16',wfm_SAR-2).',16).';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        case {'SIR_FBR_SARIN'}
            
            
            %%%%%%%%% SIN FBR Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+(i-1).*record_size,'bof');
            CS.FBR.data(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*int8=>int8'],4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+n_points+(i-1).*record_size,'bof');
            CS.FBR.N_pulses(:,i)=fread(fid,N_block,'uint16',wfm_SIN-2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+n_points+2+(i-1).*record_size,'bof');
            CS.FBR.FLAG(:,:,i)=dec2bin(fread(fid,N_block,'uint16',wfm_SIN-2).',16).';
        
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
        case {'SIR_L1B_SARIN'}
                        
            
            %%%%%%%%% SIN Average Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*0+(i-1).*record_size,'bof');
            CS.AVG.TAI.days(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*1+(i-1).*record_size,'bof');
            CS.AVG.TAI.secs(:,i)=fread(fid,1,'uint32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*2+(i-1).*record_size,'bof');
            CS.AVG.TAI.microsecs(:,i)=fread(fid,1,'uint32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*3+(i-1).*record_size,'bof');
            CS.AVG.lat(:,i)=fread(fid,1,'int32').*1e-7;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*4+(i-1).*record_size,'bof');
            CS.AVG.lon(:,i)=fread(fid,1,'int32').*1e-7;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*5+(i-1).*record_size,'bof');
            CS.AVG.H(:,i)=fread(fid,1,'int32')*1e-3;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+(i-1).*record_size,'bof');
            CS.AVG.win_delay(:,i)=fread(fid,1,'uint64')*1e-12;
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+(i-1).*record_size,'bof');
            CS.AVG.data(:,i)=fread(fid,n_points,'uint16');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+2.*n_points+(i-1).*record_size,'bof');
            CS.AVG.echo_scaling(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4.*6+8+2.*n_points+4.*1+(i-1).*record_size,'bof');
            CS.AVG.echo_scale_power(:,i)=fread(fid,1,'int32');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4*6+8+2.*n_points+4.*2+(i-1).*record_size,'bof');
            CS.AVG.N_averaged_echoes(:,i)=fread(fid,1,'uint16');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+4*6+8+2.*n_points+4.*2+2.*1+(i-1).*record_size,'bof');
            dummy=dec2bin(fread(fid,1,'uint16').',16);
            CS.AVG.OneHz_Echo_Err(:,i)=bin2dec(dummy(1).');
            CS.AVG.Mispointing_Err(:,i)=bin2dec(dummy(16).');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%% SIN Full Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+(i-1).*record_size,'bof');
            CS.SIN.data(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*uint16'],4.*2+2.*2+50*2+N_samples.*4*2+N_samples.*4*4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+(i-1).*record_size,'bof');
            CS.SIN.echo_scaling(:,i)=fread(fid,N_block,'int32',wfm_SIN-4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+4.*1+(i-1).*record_size,'bof');
            CS.SIN.echo_scale_power(:,i)=fread(fid,N_block,'int32',wfm_SIN-4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+4*2+(i-1).*record_size,'bof');
            CS.SIN.N_averaged_echoes(:,i)=fread(fid,N_block,'int16',wfm_SIN-2);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+4.*2+2.*1+(i-1).*record_size,'bof');
            CS.SIN.FLAG(:,:,i)=dec2bin(fread(fid,N_block,'uint16',wfm_SIN-2).',16).';
            CS.SIN.FLAG_tab.Approximate_Beam_Steering(:,i)=bin2dec(CS.SIN.FLAG(1,:,i).');
            CS.SIN.FLAG_tab.Exact_Beam_Steering(:,i)=bin2dec(CS.SIN.FLAG(2,:,i).');
            CS.SIN.FLAG_tab.Doppler_Weighting_Computed(:,i)=bin2dec(CS.SIN.FLAG(3,:,i).');
            CS.SIN.FLAG_tab.Doppler_Weighting_Applied_Before_Stack(:,i)=bin2dec(CS.SIN.FLAG(4,:,i).');
            CS.SIN.FLAG_tab.Multilook_Incomplete(:,i)=bin2dec(CS.SIN.FLAG(5,:,i).');
            CS.SIN.FLAG_tab.Beam_Angle_Steering_Err(:,i)=bin2dec(CS.SIN.FLAG(6,:,i).');
            CS.SIN.FLAG_tab.AntiAliased_Power_Echo(:,i)=bin2dec(CS.SIN.FLAG(7,:,i).');
            CS.SIN.FLAG_tab.Auto_Beam_Steering(:,i)=bin2dec(CS.SIN.FLAG(8,:,i).');
            CS.SIN.FLAG_tab.Reserved(:,i)=bin2dec(CS.SIN.FLAG(9:16,:,i).');
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+4*2+2.*2+(i-1).*record_size,'bof');
            CS.SIN.beam_param(:,i)=fread(fid,50.*N_block,[num2str(50) '*uint16'],wfm_SIN-2.*50);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+4*2+2.*2+50*2+(i-1).*record_size,'bof');
            CS.SIN.coherence(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*uint16'],wfm_SIN-2.*N_samples.*4);
            fseek(fid,MPH_size+HDR.SPH_SIZE+offset+av_wfm_SIN+2.*n_points+4*2+2.*2+50*2+N_samples.*4*2+(i-1).*record_size,'bof');
            CS.SIN.phase_difference(:,i)=fread(fid,n_points.*N_block,[num2str(n_points) '*int32'],wfm_SIN-4.*N_samples.*4);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    end
    
end

CS.GEO.V.V=sqrt(CS.GEO.V.Vx.^2+CS.GEO.V.Vy.^2+CS.GEO.V.Vz.^2);
CS.GEO.Serial_Sec_Num=CS.GEO.TAI.days.*24.*60.*60+CS.GEO.TAI.secs+CS.GEO.TAI.microsecs./1e6;

CS.COR.TOTAL_gim=CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.inv_bar+CS.COR.dac+CS.COR.gim_ion+CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+...
CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide;

CS.COR.TOTAL_model=CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.inv_bar+CS.COR.dac+CS.COR.model_ion+CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+...
CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide;

CS.GEO.Elapsed_Time=zeros(size(CS.GEO.Serial_Sec_Num));
CS.GEO.Elapsed_Time(CS.GEO.Serial_Sec_Num~=0)=CS.GEO.Serial_Sec_Num(CS.GEO.Serial_Sec_Num~=0)-CS.GEO.Start_Time;

CS.GEO.MODE_ID=CS.GEO.MODE_ID_Tab;
CS.GEO=rmfield(CS.GEO,'MODE_ID_Tab');

CS.GEO.INS_CFG=CS.GEO.INS_CFG_Tab;
CS.GEO=rmfield(CS.GEO,'INS_CFG_Tab');

CS.GEO.MCD_FLAG=CS.GEO.MCD_FLAG_Tab;
CS.GEO=rmfield(CS.GEO,'MCD_FLAG_Tab');

CS.COR.corr_status=CS.COR.corr_status_tab;
CS.COR=rmfield(CS.COR,'corr_status_tab');

CS.COR.corr_error=CS.COR.corr_error_tab;
CS.COR=rmfield(CS.COR,'corr_error_tab');

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'} 
        
        CS.LRM.data=reshape(CS.LRM.data,N_samples,N_block,n_recs);
        
    case {'SIR_L1B_SAR'}
        
        CS.AVG.Serial_Sec_Num=CS.AVG.TAI.days.*24.*60.*60+CS.AVG.TAI.secs+CS.AVG.TAI.microsecs./1e6;
        CS.AVG.Elapsed_Time=CS.AVG.Serial_Sec_Num-CS.GEO.Start_Time;
        CS.SAR.data=reshape(CS.SAR.data,N_samples,N_block,n_recs);
        CS.SAR.beam_param=reshape(CS.SAR.beam_param,50,N_block,n_recs);
        CS.SAR.beam_param(1,:,:)=CS.SAR.beam_param(1,:,:)./100;
        CS.SAR.beam_param(2,:,:)=CS.SAR.beam_param(2,:,:)./100;
        CS.SAR.beam_param(3,:,:)=CS.SAR.beam_param(3,:,:);
        CS.SAR.beam_param(4,:,:)=CS.SAR.beam_param(4,:,:)./100;
        CS.SAR.beam_param(5,:,:)=CS.SAR.beam_param(5,:,:).*100;
        CS.SAR.beam_param(6,:,:)=CS.SAR.beam_param(6,:,:).*1e-6;
        CS.SAR.beam_param(7,:,:)=CS.SAR.beam_param(7,:,:).*1e-6;
        CS.SAR.FLAG=CS.SAR.FLAG_tab;
        CS.SAR=rmfield(CS.SAR,'FLAG_tab');
        
    case {'SIR_L1B_SARIN'}
        
        CS.AVG.Serial_Sec_Num=CS.AVG.TAI.days.*24.*60.*60+CS.AVG.TAI.secs+CS.AVG.TAI.microsecs./1e6;
        CS.AVG.Elapsed_Time=CS.AVG.Serial_Sec_Num-CS.GEO.Start_Time;
        CS.SIN.data=reshape(CS.SIN.data,n_points,N_block,n_recs);
        CS.SIN.coherence=reshape(CS.SIN.coherence,n_points,N_block,n_recs);
        CS.SIN.phase_difference=reshape(CS.SIN.phase_difference,n_points,N_block,n_recs);
        CS.SIN.data=reshape(CS.SIN.data,n_points,N_block,n_recs);
        CS.SIN.beam_param=reshape(CS.SIN.beam_param,50,N_block,n_recs);
        CS.SIN.beam_param(1,:,:)=CS.SIN.beam_param(1,:,:)./100;
        CS.SIN.beam_param(2,:,:)=CS.SIN.beam_param(2,:,:)./100;
        CS.SIN.beam_param(3,:,:)=CS.SIN.beam_param(3,:,:);
        CS.SIN.beam_param(4,:,:)=CS.SIN.beam_param(4,:,:)./100;
        CS.SIN.beam_param(5,:,:)=CS.SIN.beam_param(5,:,:)./100;
        CS.SIN.beam_param(6,:,:)=CS.SIN.beam_param(6,:,:).*1e-6;
        CS.SIN.beam_param(7,:,:)=CS.SIN.beam_param(7,:,:).*1e-6;
        CS.SIN.FLAG=CS.SIN.FLAG_tab;
        CS.SIN=rmfield(CS.SIN,'FLAG_tab');
        
    case  {'SIR_FBR_SAR'}
        
        CS.FBR.I=CS.FBR.data(1:2:end,:);
        CS.FBR.Q=CS.FBR.data(2:2:end,:);
        CS.FBR.data=reshape((complex(CS.FBR.I,CS.FBR.Q)),N_samples,[]);
        CS.FBR=rmfield(CS.FBR, {'I';'Q'});
        
    case {'SIR_FBR_SARIN'}
        
        I = CS.FBR.data(1:2:end,:);
        Q = CS.FBR.data(2:2:end,:);
        Ib = reshape(I,N_SIN_samples,SAR_pulses_burst,2,N_block,n_recs);
        Qb = reshape(Q,N_SIN_samples,SAR_pulses_burst,2,N_block,n_recs);
                
        CS.FBR.Rx1 = squeeze(complex(Ib(:,:,1,:,:),Qb(:,:,1,:,:)));
        CS.FBR.Rx2 = squeeze(complex(Ib(:,:,2,:,:),Qb(:,:,2,:,:)));
        
        clear tmp Rx1 Rx2;
        CS.FBR=rmfield(CS.FBR,{'data'});
end


fclose(fid);


