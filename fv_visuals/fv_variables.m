% /////// fv_variables ///////
% Simple function with 2 uses
% 1: Generate a list of all TUFLOW-FV variable names (case sensitive)
%    list = fv_variables
% 2: input a list of variable names and output a cell array of variable names which can be found in a TUFLOW-FV netcdf output file
%    If you specify a variable which is unrecogniseable then an error is thrown avoiding confusing errors down the track
%    varnames_out = fv_variables(varnames_in)
%
% example (#2):
%   varnames_out = fv_variables({'H';'V_mag';'tss'});
%   varnames_out = {'H';'V_x';'V_y';'TSS'};
%
% JN July 2012

function varnames_out = fv_variables(varargin)

% are we generating a list of all variables
if isempty(varargin)
    varnames_out = fv_variables_list;
    return
else
    varnames_in = varargin{:};
    if ~iscell(varnames_in)
        if size(varnames_in,1) == 1
            varnames_in = {varnames_in};
        else
            error('expecting cell array for input varnames_in')
        end
    end
end


% do variables need to be split into components
nvi = length(varnames_in);
skip_v=false;
skip_w=false;
k=1;
for aa = 1:nvi
    switch lower(varnames_in{aa})
        case {'v','v_mag','vmag','v_dir','vdir','current'}
            if ~skip_v
                varnames_out{k,1} = 'V_x';
                varnames_out{k+1,1} = 'V_y';
                k = k+2;
                skip_v = true;
            else continue
            end
        case {'w10','w10_mag','w10mag','w10_dir','w10dir'}
            if ~skip_w
                varnames_out{k,1} = 'W10_x';
                varnames_out{k+1,1} = 'W10_y';
                k = k+2;
                skip_w = true;
            else continue
            end
        case {'wvstr_mag','wvstr_dir'}
            varnames_out{k,1} = 'WVSTR_x';
            varnames_out{k+1,1} = 'WVSTR_y';
            k = k+2;
        case {'sedload_total'}
            varnames_out{k,1} = 'SEDLOAD_TOTAL_x';
            varnames_out{k+1,1} = 'SEDLOAD_TOTAL_y';
            k = k+2;
        case {'turb'}
            varnames_out{k,1} = 'TSS';
            k=k+1;
        case {'bedload_total'}
            varnames_out{k,1} = 'BEDLOAD_TOTAL_x';
            varnames_out{k+1,1} = 'BEDLOAD_TOTAL_y';
            k = k+2;
        case {'suspload_total'}
            varnames_out{k,1} = 'SUSPLOAD_TOTAL_x';
            varnames_out{k+1,1} = 'SUSPLOAD_TOTAL_y';
            k = k+2;
        otherwise
            varnames_out{k,1}=varnames_in{aa};
            k=k+1;
    end
end
nvo = length(varnames_out);

% are variables names or similar name (same leters but different case) legal TUFLOW-FV variable names
list = fv_variables_list;
for aa = 1:nvo
    v_name = varnames_out{aa};
    if ismember(v_name,list)
        continue
    else
        [loj i] = ismember(lower(v_name),lower(list));
        if loj
            v_name = list{i};
            varnames_out{aa} = v_name;
            continue
        else
            i = strfind(v_name,'_');
            if ~isempty(i)
                i = i(end);
                numstr = v_name(i+1:end);
                num = str2double(numstr);
                if ~isnan(num)
                    v_name = strrep(v_name,numstr,'#');
                    [loj i] = ismember(lower(v_name),lower(list));
                    if loj
                        v_name = list{i};
                        v_name = strrep(v_name,'#',numstr);
                        varnames_out{aa} = v_name;
                        continue
                    end
                end
            end
        end
        error([v_name ' is not a recognised TUFLOW-FV variable name'])
    end
end

varnames_out = unique(varnames_out);

% /////// nested functions ///////
% list = fv_variables
%
% When called fv_variables provides a cell array of all the variable names
% known to TUFLOW-FV.

function list = fv_variables_list

hydro = {'H';'D';'V_x';'V_y';'W'};

wind = {'W10_x';'W10_y'};

atmo = {'MSLP';'AIR_TEMP';'REL_HUM';'SW_RAD';'LW_RAD';'PRECIP';'EVAP'};

wave = {'WVHT';'WVPER';'WVDIR';'WVSTR_x';'WVSTR_y';'WVUBOT';'WVPERBOT'};

scalar =  {'TEMP';'SAL';'PTM_#'};

sedi = {'TSS';'SED_#';'TURB'};

trace = {'TRACE_#';'PTM_#'};

bathy = {'ZB';'DZB';'RHOW'};

turbz = {'TURBZ_TKE';'TURBZ_EPS';'TURBZ_L';'TURBZ_SPFSQ';'TURBZ_BVFSQ';...
    'TURBZ_NUM';'TURBZ_NUH';'TURBZ_NUS';};

bed = {'TAU';'TAUB';'TAUC';'TAUW';'TAUCW';'PICKUP';'PICKUP_TOTAL';...
    'DEPOSITION';'DEPOSITION_TOTAL';'NETSEDRATE_TOTAL';'BED_MASS';'BED_MASS_TOTAL';...
    'BED_MASS_SED_#';'BED_MASS_LAYER_#';...
    'BEDLOAD_TOTAL_x';'BEDLOAD_TOTAL_y';...
    'SUSPLOAD_TOTAL_x';'SUSPLOAD_TOTAL_y';...
    'SEDLOAD_TOTAL_x';'SEDLOAD_TOTAL_y';'THICK'};

aed = {'WQ_AED_NITROGEN_AMM';'WQ_AED_NITROGEN_NIT';...
    'WQ_AED_ORGANIC_MATTER_DOC';'WQ_AED_ORGANIC_MATTER_DON';...
    'WQ_AED_ORGANIC_MATTER_DOP';'WQ_AED_ORGANIC_MATTER_POC';...
    'WQ_AED_ORGANIC_MATTER_PON';'WQ_AED_ORGANIC_MATTER_POP';...
    'WQ_AED_OXYGEN_OXY';'WQ_AED_PHOSPHORUS_FRP';...
    'WQ_AED_PHOSPHORUS_FRP_ADS';
    'WQ_AED_PHYTOPLANKTON_GRN';...
    'WQ_AED_PHYTOPLANKTON_GRN_IN';'WQ_AED_PHYTOPLANKTON_GRN_IP';...
    'WQ_AED_PHYTOPLANKTON_BGA';...
    'WQ_AED_PHYTOPLANKTON_BGA_IN';'WQ_AED_PHYTOPLANKTON_BGA_IP';...
    'WQ_AED_PHYTOPLANKTON_FDIAT';...
    'WQ_AED_PHYTOPLANKTON_FDIAT_IN';'WQ_AED_PHYTOPLANKTON_FDIAT_IP';...
    'WQ_AED_PHYTOPLANKTON_MDIAT';...
    'WQ_AED_PHYTOPLANKTON_MDIAT_IN';'WQ_AED_PHYTOPLANKTON_MDIAT_IP';...
    'WQ_AED_PHYTOPLANKTON_CRYPTO';...
    'WQ_AED_PHYTOPLANKTON_CRYPTO_IN';'WQ_AED_PHYTOPLANKTON_CRYPTO_IP';...
    'WQ_AED_PHYTOPLANKTON_DIATOM';'WQ_AED_PHYTOPLANKTON_DIATOM_IN';...
    'WQ_AED_PHYTOPLANKTON_DIATOM_IP';'WQ_AED_PHYTOPLANKTON_GREEN';...
    'WQ_AED_PHYTOPLANKTON_GREEN_IN';'WQ_AED_PHYTOPLANKTON_GREEN_IP';...
    'WQ_AED_PHYTOPLANKTON_PHY01';'WQ_AED_PHYTOPLANKTON_PHY02';...
    'WQ_AED_PHYTOPLANKTON_PHY03';'WQ_AED_PHYTOPLANKTON_PHY04';...
    'WQ_AED_PHYTOPLANKTON_PHY05';'WQ_AED_PHYTOPLANKTON_PHY06';...
    'WQ_AED_PHYTOPLANKTON_PHY07';'WQ_AED_PHYTOPLANKTON_PHY08';...
    'WQ_AED_PHYTOPLANKTON_PHY09';'WQ_AED_PHYTOPLANKTON_PHY10';...
    'WQ_AED_PHYTOPLANKTON_PHY11';'WQ_AED_SILICA_RSI';...
    'WQ_DIAG_AED_PHYTOPLANKTON_GRN_FNIT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_GRN_FPHO';...
    'WQ_DIAG_AED_PHYTOPLANKTON_GRN_FSIL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_GRN_FT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_GRN_FI';...
    'WQ_DIAG_AED_PHYTOPLANKTON_GRN_FSAL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_BGA_FNIT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_BGA_FPHO';...
    'WQ_DIAG_AED_PHYTOPLANKTON_BGA_FSIL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_BGA_FT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_BGA_FI';...
    'WQ_DIAG_AED_PHYTOPLANKTON_BGA_FSAL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_FDIAT_FNIT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_FDIAT_FPHO';...
    'WQ_DIAG_AED_PHYTOPLANKTON_FDIAT_FSIL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_FDIAT_FT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_FDIAT_FI';...
    'WQ_DIAG_AED_PHYTOPLANKTON_FDIAT_FSAL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_MDIAT_FNIT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_MDIAT_FPHO';...
    'WQ_DIAG_AED_PHYTOPLANKTON_MDIAT_FSIL';...
    'WQ_DIAG_AED_PHYTOPLANKTON_MDIAT_FT';...
    'WQ_DIAG_AED_PHYTOPLANKTON_MDIAT_FI';...
    'WQ_DIAG_AED_PHYTOPLANKTON_MDIAT_FSAL';...
    'WQ_DIAG_AED_NITROGEN_DENIT';'WQ_DIAG_AED_NITROGEN_NITRIF';...
    'WQ_DIAG_AED_NITROGEN_SED_AMM';'WQ_AED_PATHOGENS_CRYPT';'WQ_AED_PATHOGENS_ECOLI';...
    'WQ_AED_PATHOGENS_FCOLI';'WQ_AED_PATHOGENS_ENT';...
    'WQ_DIAG_AED_NITROGEN_SED_NIT';...
    'WQ_DIAG_AED_ORGANIC_MATTER_BOD';'WQ_DIAG_AED_ORGANIC_MATTER_DOC_MINER';...
    'WQ_DIAG_AED_ORGANIC_MATTER_DON_MINER';'WQ_DIAG_AED_ORGANIC_MATTER_DOP_MINER';...
    'WQ_DIAG_AED_ORGANIC_MATTER_POC_MINER';'WQ_DIAG_AED_ORGANIC_MATTER_PON_MINER';...
    'WQ_DIAG_AED_ORGANIC_MATTER_POP_MINER';'WQ_DIAG_AED_ORGANIC_MATTER_SED_DOC';...
    'WQ_DIAG_AED_ORGANIC_MATTER_SED_DON';'WQ_DIAG_AED_ORGANIC_MATTER_SED_DOP';...
    'WQ_DIAG_AED_ORGANIC_MATTER_SED_POC';'WQ_DIAG_AED_ORGANIC_MATTER_SED_PON';...
    'WQ_DIAG_AED_ORGANIC_MATTER_SED_POP';'WQ_DIAG_AED_OXYGEN_AED_OXYGEN_SAT';...
    'WQ_DIAG_AED_OXYGEN_ATM_OXY_EXCH';'WQ_DIAG_AED_OXYGEN_SED_OXY';...
    'WQ_DIAG_AED_PHOSPHORUS_SED_FRP';'WQ_DIAG_AED_PHYTOPLANKTON_GPP';...
    'WQ_DIAG_AED_PHYTOPLANKTON_IN';'WQ_DIAG_AED_PHYTOPLANKTON_IP';...
    'WQ_DIAG_AED_PHYTOPLANKTON_NCP';'WQ_DIAG_AED_PHYTOPLANKTON_NPR';...
    'WQ_DIAG_AED_PHYTOPLANKTON_PAR';'WQ_DIAG_AED_PHYTOPLANKTON_PPR';...
    'WQ_DIAG_AED_PHYTOPLANKTON_TCHLA';'WQ_DIAG_AED_SILICA_SED_RSI';...
    'WQ_DIAG_AED_TOTALS_AED_TOTALS_TN';'WQ_DIAG_AED_TOTALS_AED_TOTALS_TOC';...
    'WQ_DIAG_AED_TOTALS_AED_TOTALS_TP';'WQ_DIAG_AED_TOTALS_AED_TOTALS_TSS';...
    'WQ_DIAG_AED_TOTALS_AED_TOTALS_TURBIDIT';'WQ_OXY_OXY';'WQ_CAR_RSI';'WQ_NIT_AMM';...
    'WQ_NIT_NIT';'WQ_PHS_FRP';'WQ_PHS_FRP_ADS';'WQ_OGM_DON';'WQ_OGM_PON';'WQ_OGM_DOP';...
    'WQ_OGM_POP';'WQ_OGM_DOC';'WQ_OGM_POC';'WQ_PHY_GRN';'WQ_DIAG_OXY_SED_OXY';...
    'WQ_DIAG_OXY_ATM_OXY_EXCH';'WQ_DIAG_OXY_AED_OXYGEN_SAT';'WQ_DIAG_CAR_SED_RSI';...
    'WQ_DIAG_NIT_NITRIF';'WQ_DIAG_NIT_DENIT';'WQ_DIAG_NIT_SED_AMM';'WQ_DIAG_NIT_SED_NIT';...
    'WQ_DIAG_PHS_SED_FRP';'WQ_DIAG_OGM_PON_MINER';'WQ_DIAG_OGM_DON_MINER';...
    'WQ_DIAG_OGM_SED_PON';'WQ_DIAG_OGM_SED_DON';'WQ_DIAG_OGM_POP_MINER';...
    'WQ_DIAG_OGM_DOP_MINER';'WQ_DIAG_OGM_SED_POP';'WQ_DIAG_OGM_SED_DOP';...
    'WQ_DIAG_OGM_POC_MINER';'WQ_DIAG_OGM_DOC_MINER';'WQ_DIAG_OGM_SED_POC';...
    'WQ_DIAG_OGM_SED_DOC';'WQ_DIAG_OGM_BOD';'WQ_DIAG_PHY_GRN_NTOP';'WQ_DIAG_PHY_GRN_FT';...
    'WQ_DIAG_PHY_GRN_FI';'WQ_DIAG_PHY_GRN_FNIT';'WQ_DIAG_PHY_GRN_FPHO';...
    'WQ_DIAG_PHY_GRN_FSIL';'WQ_DIAG_PHY_GRN_FSAL';'WQ_DIAG_PHY_GPP';'WQ_DIAG_PHY_NCP';...
    'WQ_DIAG_PHY_PPR';'WQ_DIAG_PHY_NPR';'WQ_DIAG_PHY_NUP';'WQ_DIAG_PHY_PUP';...
    'WQ_DIAG_PHY_CUP';'WQ_DIAG_PHY_PAR';'WQ_DIAG_PHY_TCHLA';'WQ_DIAG_PHY_TPHYS';...
    'WQ_DIAG_PHY_IN';'WQ_DIAG_PHY_IP';'WQ_DIAG_TOT_AED_TOTALS_TN';'WQ_DIAG_TOT_AED_TOTALS_TP';...
    'WQ_DIAG_TOT_AED_TOTALS_TOC';'WQ_DIAG_TOT_AED_TOTALS_TSS';'WQ_DIAG_TOT_AED_TOTALS_TURBIDITY';...
    'WQ_PHY_MDIAT';'WQ_PHY_MDIAT';'WQ_PHY_FDIAT';'WQ_PHY_BGA';'WQ_PHY_CLD';
    'WQ_PAT_ENT';'WQ_PAT_FCOLI';'WQ_PAT_EC';...
    'WQ_OXY_OXY';'WQ_CAR_RSI';'WQ_NIT_AMM';'WQ_NIT_NIT';...
    'WQ_PHS_FRP';'WQ_PHS_FRP_ADS';'WQ_OGM_DON';...
    'WQ_OGM_PON';'WQ_OGM_DOP';'WQ_OGM_POP';...
    'WQ_OGM_DOC';'WQ_OGM_POC';'WQ_PHY_GRN';...
    'WQ_PHY_BGA';'WQ_PHY_FDIAT';'WQ_PHY_MDIAT';...
    'WQ_PAT_ECOLI';'WQ_PAT_FCOLI';'WQ_PAT_ENT';...
    'WQ_TRC_SS1';'WQ_SIL_RSI';...
    'WQ_DIAG_SDF_FSED_OXY';'WQ_DIAG_OXY_SED_OXY';...
    'WQ_DIAG_OXY_ATM_OXY_EXCH';...
    'WQ_DIAG_OXY_SAT';'WQ_DIAG_OGM_BOD';...
    'WQ_DIAG_NIT_NITRIF';'WQ_DIAG_OGM_DOC_MINER';...
    'WQ_DIAG_PHY_CUP';'WQ_DIAG_SDF_FSED_RSI';'WQ_DIAG_SDF_FSED_AMM';...
    'WQ_DIAG_SDF_FSED_NIT';'WQ_DIAG_SDF_FSED_FRP';'WQ_DIAG_SDF_FSED_PON';...
    'WQ_DIAG_SDF_FSED_DON';'WQ_DIAG_SDF_FSED_POP';'WQ_DIAG_SDF_FSED_DOP';...
    'WQ_DIAG_SDF_FSED_POC';'WQ_DIAG_SDF_FSED_DOC';'WQ_DIAG_SDF_FSED_DIC';...
    'WQ_DIAG_SIL_SED_RSI';'WQ_DIAG_NIT_DENIT';'WQ_DIAG_NIT_SED_AMM';...
    'WQ_DIAG_NIT_SED_NIT';'WQ_DIAG_PHS_SED_FRP';'WQ_DIAG_OGM_PON_MINER';...
    'WQ_DIAG_OGM_DON_MINER';'WQ_DIAG_OGM_SED_PON';'WQ_DIAG_OGM_SED_DON';...
    'WQ_DIAG_OGM_POP_MINER';'WQ_DIAG_OGM_DOP_MINER';'WQ_DIAG_OGM_SED_POP';...
    'WQ_DIAG_OGM_SED_DOP';'WQ_DIAG_OGM_POC_MINER';'WQ_DIAG_OGM_SED_POC';...
    'WQ_DIAG_OGM_SED_DOC';'WQ_DIAG_PHY_GRN_NTOP';'WQ_DIAG_PHY_GPP';...
    'WQ_DIAG_PHY_NCP';'WQ_DIAG_PHY_PPR';'WQ_DIAG_PHY_NPR';'WQ_DIAG_PHY_NUP';...
    'WQ_DIAG_PHY_PUP';'WQ_DIAG_PHY_PAR';'WQ_DIAG_PHY_TCHLA';...
    'WQ_DIAG_PHY_TPHYS';'WQ_DIAG_PHY_IN';'WQ_DIAG_PHY_IP';...
    'WQ_DIAG_TOT_TN';'WQ_DIAG_TOT_TP';'WQ_DIAG_TOT_TOC';'WQ_DIAG_TOT_TSS';'WQ_DIAG_TOT_TURBIDITY';
    'WQ_PAT_CRYPTO';...
		'WQ_PAT_ECOLI';...
		'WQ_PAT_ENTERO';...
		'WQ_PAT_CAMPY';...
		'WQ_PAT_VIRUS';... 
    'PAR';
    };

list = cat(1,hydro,wind,atmo,wave,scalar,sedi,trace,bathy,turbz,bed,aed);