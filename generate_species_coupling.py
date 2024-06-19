############################################
# Generate species coupling code.
# Used for GEOS-Chem with CESM and WRF-GC
# (c) 2022-2024 Haipeng Lin <myself@jimmielin.me>
#
# Licensed under the GNU General Public License v2
############################################

import re
from operator import itemgetter

# User configurable options:

# choose input file:
input_file_name = "v14.4.0.yml"

# cesm | wrfgc
mode = "cesm"

# choose parser_mode: "regex" for 13.4.x, or "yaml" for 14.1.0+
parser_mode = "yaml" if "yml" in input_file_name else "regex"

# is version 14 and above? used for WRF-GC
# if yes, uses the new format (State_Chm%Species(id_Hg2)%Conc(:,:,:)) for State_Chm%
is_gc_14_and_above = "14" in input_file_name

#
# These are "quirky" species, which necessitate special handling
# in coupled models such as WRF and CESM.
#

# Non-gas-phase species ("Aerosols")
nongas = ["AERI", "ASOAN", "ASOA1", "ASOA2", "ASOA3", "ASOG1", "ASOG2", "ASOG3", "AONITA", "BCPI", "BCPO", "BrSALA", "BrSALC", "DMS", "DST1", "DST2", "DST3", "DST4", "INDIOL", "IONITA", "ISALA", "ISALC", "ISN1OA", "ISN1OG", "ISOA1", "ISOA2", "ISOA3", "MONITA", "MSA", "NH4", "NIT", "NITs", "OCPI", "OCPO", "OPOA1", "OPOA2", "OPOG1", "OPOG2", "SALA", "SALC", "SALACl", "SALAAL", "SALCAL", "SALCCl", "SO4", "SO4s", "SOAIE", "SOAGX", "SOAME", "SOAMG", "SOAP", "SOAS", "TSOA0", "TSOA1", "TSOA2", "TSOA3", "pFe", "POA1", "POA2", "POG1", "POG2"]
handled_by_mam4 = ["SO4", "OCPI", "OCPO", "BCPI", "BCPO", "DST1", "DST2", "DST3", "DST4", "SALA", "SALC"]

# ComplexSOA
complexsoa = ["ASOA1", "ASOA2", "ASOA3", "ASOAN", "ASOG1", "ASOG2", "ASOG3", "TSOA0", "TSOA1", "TSOA2", "TSOA3", "TSOG0", "TSOG1", "TSOG2", "TSOG3"]
complexsoa_svpoa = ["NAP", "POA1", "POA2", "POG1", "POG2", "OPOA1", "OPOA2", "OPOG1", "OPOG2"]
simplesoa = ["SOAS", "SOAP"]

###################
#  --- CESM ---   #
###################
# CESM has a few fixed species that need to be written to inv(ariants) list
# from gckpp_Parameters.F90, this is H2, N2, O2, RCOOH
cesm_invariants = ["H2", "N2", "O2", "RCOOH"]

# CESM-GC requires extra species for aerosols
cesm_aer = ['bc_a1','bc_a4','dst_a1','dst_a2','dst_a3','ncl_a1','ncl_a2','ncl_a3','num_a1','num_a2','num_a3','num_a4','pom_a1','pom_a4','so4_a1','so4_a2','so4_a3','soa1_a1','soa1_a2','soa2_a1','soa2_a2','soa3_a1','soa3_a2','soa4_a1','soa4_a2','soa5_a1','soa5_a2','H2SO4','SOAG0','SOAG1','SOAG2','SOAG3','SOAG4']
cesm_aer_mass = [12.011000, 12.011000, 135.064039, 135.064039, 135.064039, 58.442468, 58.442468, 58.442468, 1.007400, 1.007400, 1.007400, 1.007400, 12.011000, 12.011000, 115.107340, 115.107340, 115.107340, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000, 98.078400, 250.445000, 250.445000, 250.445000, 250.445000, 250.445000]

###################
# --- WRF-GC ---  #
###################
# WRF requires non-gas species to be separated into the end. use nongas.

# WRF-GC v2.0+ (Feng et al., 2021) also requires diagnostic bins for aerosols
wrf_extra_coupling = "diag_so4_a1,diag_so4_a2,diag_so4_a3,diag_so4_a4,diag_nit_a1,diag_nit_a2,diag_nit_a3,diag_nit_a4,diag_nh4_a1,diag_nh4_a2,diag_nh4_a3,diag_nh4_a4,diag_ocpi_a1,diag_ocpi_a2,diag_ocpi_a3,diag_ocpi_a4,diag_ocpo_a1,diag_ocpo_a2,diag_ocpo_a3,diag_ocpo_a4,diag_bcpi_a1,diag_bcpi_a2,diag_bcpi_a3,diag_bcpi_a4,diag_bcpo_a1,diag_bcpo_a2,diag_bcpo_a3,diag_bcpo_a4,diag_seas_a1,diag_seas_a2,diag_seas_a3,diag_seas_a4,diag_dst_a1,diag_dst_a2,diag_dst_a3,diag_dst_a4,diag_soas_a1,diag_soas_a2,diag_soas_a3,diag_soas_a4,diag_so4_cw1,diag_so4_cw2,diag_so4_cw3,diag_so4_cw4,diag_nit_cw1,diag_nit_cw2,diag_nit_cw3,diag_nit_cw4,diag_nh4_cw1,diag_nh4_cw2,diag_nh4_cw3,diag_nh4_cw4,diag_ocpi_cw1,diag_ocpi_cw2,diag_ocpi_cw3,diag_ocpi_cw4,diag_ocpo_cw1,diag_ocpo_cw2,diag_ocpo_cw3,diag_ocpo_cw4,diag_bcpi_cw1,diag_bcpi_cw2,diag_bcpi_cw3,diag_bcpi_cw4,diag_bcpo_cw1,diag_bcpo_cw2,diag_bcpo_cw3,diag_bcpo_cw4,diag_seas_cw1,diag_seas_cw2,diag_seas_cw3,diag_seas_cw4,diag_dst_cw1,diag_dst_cw2,diag_dst_cw3,diag_dst_cw4,diag_soas_cw1,diag_soas_cw2,diag_soas_cw3,diag_soas_cw4,diag_water_a1,diag_water_a2,diag_water_a3,diag_water_a4,diag_num_a1,diag_num_a2,diag_num_a3,diag_num_a4,diag_num_cw1,diag_num_cw2,diag_num_cw3,diag_num_cw4"

# Starting in WRF-GC v3.0 (WRF v4 + GEOS-Chem v13.4.1/v14.0.0), GEOS-Chem species
# names will directly be used in WRF after a lowercase processing. This is to
# resolve many lingering headaches. (hplin, 8/11/22)

############ NO USER CONFIGURABLE CODE BELOW ############

if parser_mode == "regex":
    with open(input_file_name, 'r') as file:
        input_gc = file.read()
        
    out_regex = r"N:\s*(\d+)\s*(\w+)\s*Adv:\s*(F|T)\sDd:\s*(F|T)\sWd:\s*(F|T)\s*MW_g:\s*(\d+\.\d+)"
    parsed_spc_list = re.findall(out_regex, input_gc)
elif parser_mode == "yaml":
    import yaml

    with open(input_file_name, 'r') as file:
        try:
            input_gc_obj = yaml.safe_load(file)
            parsed_spc_list = []
            N = 0

            # dump it into 5-tuples. no need to sort by advect - it will be sorted later on
            for spc in input_gc_obj:
                # KLUDGE: for some reason pyyaml parses 'NO' == False ???
                if not spc:
                    spcName = 'NO'
                else:
                    spcName = spc

                N += 1
                Is_Advected = 'T' if input_gc_obj[spc].get("Is_Advected", False) else 'F'
                Is_DryDep = 'T' if input_gc_obj[spc].get("Is_DryDep", False) else 'F'
                Is_WetDep = 'T' if input_gc_obj[spc].get("Is_WetDep", False) else 'F'
                Mw_g = input_gc_obj[spc].get("MW_g", 0.00)


                # in non-SVPOA simulations, CESM still considers NAP a non-advect tracer
                # TODO add option if SVPOA
                if spcName == "NAP":
                    Is_Advected = 'F'

                # for certain species that GEOS-Chem reports as non-advected, still
                # report them as advected. this is because for CESM, non-advected species
                # are not "constituents" and initial conditions cannot be specified for them,
                # and they cannot be output. some special cases only here:
                # OH, HO2 because we usually want output
                # MO2, MCO3 (note GC MCO3 /= CAM-chem MCO3!!) for research purposes
                #if spcName == "OH" or spcName == "HO2":
                #    Is_Advected = 'T'

                # heap alloc.
                parsed_spc_list.append((N, spcName, Is_Advected, Is_DryDep, Is_WetDep, Mw_g))

            # print(parsed_spc_list)
        except yaml.YAMLError as exc:
            print(exc)
else:
    raise ValueError("parser_mode is neither regex nor yaml!")

# 0: N
# 1: spc name
# 2: is_adv
# 3: is_drydep
# 4: is_wetdep
# 5: Mw_g
idx_spcName = 1
idx_Advect = 2
idx_DryDep = 3
idx_WetDep = 4
idx_Mw_g = 5

# sort by adv to keep them first - even though GC should keep them first
# this is for sanity.
# parsed_spc_list.sort(key = lambda x: 0 if x[idx_Advect] == "T" else 1)

# also sort by tracer name alphabetically (but make uppercase first)
parsed_spc_list = sorted(parsed_spc_list, key = lambda x: (0 if x[idx_Advect] == "T" else 1, x[idx_spcName].upper()))

# ClOO should follow ClO, do not sort CLOCK in between
tmpA, tmpB = -1, -1
for i, x in enumerate(parsed_spc_list):
    if x[idx_spcName] == "ClOO":
        tmpA = i
    elif x[idx_spcName] == "CLOCK":
        tmpB = i
if tmpA > 0 and tmpB > 0 and tmpA > tmpB:
    parsed_spc_list[tmpA], parsed_spc_list[tmpB] = parsed_spc_list[tmpB], parsed_spc_list[tmpA]

assert len(cesm_aer) == len(cesm_aer_mass)
nongas_upper = list(map(lambda x: x.upper(), nongas))

if mode == "cesm":

    # CESM format...
    # generate solsym, adv_mass arrays
    # count total number of species first
    total_num_parsed_m1 = len(parsed_spc_list) - 1
    # print(total_num_parsed_m1)

    # first, loop through advected species which should go first
    pretty_column_counter = 0
    total_counter = 0
    nadv_counter = 0
    nonadv_counter = 0
    first_nonadv = True
    have_given_solsym_comment = False

    # Final strings. Headers will be filled later, as we need to count.
    solsym = ""
    adv_mass = ""
    drydep_list = ""
    wetdep_list = ""
    aer_drydep_list = "'dst_a1','so4_a1','nh4_a1','pom_a1','pomff1_a1','pombb1_a1','soa_a1','bc_a1','ncl_a1','num_a1','so4_a2','nh4_a2','soa_a2','ncl_a2','dst_a2','num_a2','dst_a3','ncl_a3','so4_a3','pom_a3','bc_a3','num_a3','ncl_a4','so4_a4','pom_a4','pomff1_a4','pombb1_a4','bc_a4','nh4_a4','num_a4','dst_a5','so4_a5','nh4_a5','num_a5','ncl_a6','so4_a6','nh4_a6','num_a6','dst_a7','so4_a7','nh4_a7','num_a7','soa1_a1','soa1_a2','soa2_a1','soa2_a2','soa3_a1','soa3_a2','soa4_a1','soa4_a2','soa5_a1','soa5_a2','soaff1_a1','soaff2_a1','soaff3_a1','soaff4_a1','soaff5_a1','soabb1_a1','soabb2_a1','soabb3_a1','soabb4_a1','soabb5_a1','soabg1_a1','soabg2_a1','soabg3_a1','soabg4_a1','soabg5_a1','soaff1_a2','soaff2_a2','soaff3_a2','soaff4_a2','soaff5_a2','soabb1_a2','soabb2_a2','soabb3_a2','soabb4_a2','soabb5_a2','soabg1_a2','soabg2_a2','soabg3_a2','soabg4_a2','soabg5_a2',"
    aer_wetdep_list = "'dst_a1','so4_a1','nh4_a1','pom_a1','pomff1_a1','pombb1_a1','soa_a1','bc_a1','ncl_a1','num_a1','so4_a2','nh4_a2','soa_a2','ncl_a2','dst_a2','num_a2','dst_a3','ncl_a3','so4_a3','pom_a3','bc_a3','num_a3','ncl_a4','so4_a4','pom_a4','pomff1_a4','pombb1_a4','bc_a4','nh4_a4','num_a4','dst_a5','so4_a5','nh4_a5','num_a5','ncl_a6','so4_a6','nh4_a6','num_a6','dst_a7','so4_a7','nh4_a7','num_a7','soa1_a1','soa1_a2','soa2_a1','soa2_a2','soa3_a1','soa3_a2','soa4_a1','soa4_a2','soa5_a1','soa5_a2','soaff1_a1','soaff2_a1','soaff3_a1','soaff4_a1','soaff5_a1','soabb1_a1','soabb2_a1','soabb3_a1','soabb4_a1','soabb5_a1','soabg1_a1','soabg2_a1','soabg3_a1','soabg4_a1','soabg5_a1','soaff1_a2','soaff2_a2','soaff3_a2','soaff4_a2','soaff5_a2','soabb1_a2','soabb2_a2','soabb3_a2','soabb4_a2','soabb5_a2','soabg1_a2','soabg2_a2','soabg3_a2','soabg4_a2','soabg5_a2',"

    for spc_idx, spc in enumerate(parsed_spc_list):
        # invariants do not need to be skipped, they are actually duplicated
        # if spc[1] in cesm_invariants:
        #     continue

        # CESM-GC always uses Complex SOA.
        if spc[idx_spcName] in simplesoa:
            continue

        # CESM-GC does not use semivolatile POA. For some reason, NAP is included
        # in the YML output even if it is not a SVPOA simulation.
        if spc[idx_spcName] in complexsoa_svpoa:
            # However, NAP is actually a non-advected tracer in a non-SVPOA sim
            # because it is included in KPP. So it has to be manually overridden
            if not spc[idx_spcName] == "NAP":
                continue

        # insert drydep / wetdep
        if spc[idx_DryDep] == "T":
            if (not spc[idx_spcName].upper() in nongas_upper):
                drydep_list += "'" + spc[idx_spcName].upper() + "',"
            else:
                # also, exclude species handled by GC bulk-to-MAM4 modal coupling.
                # these are overwritten at every time step and thus do not necessitate
                # any handling.
                if not spc[idx_spcName] in handled_by_mam4:
                    drydep_list += "'" + spc[idx_spcName].upper() + "'," # aer_
                    # note: temporarily, drydep of aerosols is handled in the gas list.
                    # this remains to be discussed, but is the approach taken in previous
                    # versions of CESM-GC.
        if spc[idx_WetDep] == "T":
            if (not spc[idx_spcName].upper() in nongas_upper):
                wetdep_list += "'" + spc[idx_spcName].upper() + "',"
            else:
                # also, exclude species handled by GC bulk-to-MAM4 modal coupling.
                # these are overwritten at every time step and thus do not necessitate
                # any handling.
                if not spc[idx_spcName] in handled_by_mam4:
                    wetdep_list += "'" + spc[idx_spcName].upper() + "'," # aer_
                    # note: temporarily, drydep of aerosols is handled in the gas list.
                    # this remains to be discussed, but is the approach taken in previous
                    # versions of CESM-GC.
        
        # insert aerosols in between T and F. check for first F
        if first_nonadv and spc[idx_Advect] == "F":
            first_nonadv = False
            nadv_counter = total_counter + len(cesm_aer) # ...save up to now
            # insert cesm_aer, cesm_aer_mass ...
            for aer_idx, aer_spc in enumerate(cesm_aer):
                solsym += "'" + aer_spc.ljust(15) + "', "
                adv_mass += str(cesm_aer_mass[aer_idx]).rjust(14) + "_r8, "

                pretty_column_counter += 1
                total_counter += 1
                if pretty_column_counter == 3:
                    pretty_column_counter = 0
                    if not have_given_solsym_comment:
                        solsym += "& ! Species after MAM aerosols are non-advected and will not be constituents\r\n"
                        have_given_solsym_comment = True
                    else:
                        solsym += "&\r\n"
                    adv_mass += "&\r\n"

        if spc[idx_Advect] == "F":
            nonadv_counter += 1
        
        # solsym needs to be in uppercase or FLDLST for history will complain.
        solsym += "'" + spc[idx_spcName].ljust(15).upper() + "', "
        adv_mass += str(spc[idx_Mw_g]).rjust(14) + "_r8, "
        
        pretty_column_counter += 1
        total_counter += 1
        if pretty_column_counter == 3:
            pretty_column_counter = 0
            solsym += "&\r\n"
            adv_mass += "&\r\n"

    solsym = "solsym(:" + str(total_counter) + ") = (/ &\r\n" + solsym + "/)"
    adv_mass = "adv_mass(:" + str(total_counter) + ") = (/ &\r\n" + adv_mass + "/)"

    # strip last comma - ugly code
    solsym = " ".join(solsym.rsplit(",", 1))
    adv_mass = " ".join(adv_mass.rsplit(",", 1))
    drydep_list = "".join(drydep_list.rsplit(",", 1))
    wetdep_list = "".join(wetdep_list.rsplit(",", 1))
    aer_drydep_list = "".join(aer_drydep_list.rsplit(",", 1))
    aer_wetdep_list = "".join(aer_wetdep_list.rsplit(",", 1))

    print(solsym)
    print(adv_mass)
    print("<drydep_list>" + drydep_list + "</drydep_list>")
    print("<gas_wetdep_list>" + wetdep_list + "</gas_wetdep_list>")
    print("<aer_drydep_list>" + aer_drydep_list + "</aer_drydep_list>")
    print("<aer_wetdep_list>" + aer_wetdep_list + "</aer_wetdep_list>")
    print("update chem_mods.F90: gas_pcnst = " + str(total_counter))
    print("update chem_mods.F90: nTracersMax = " + str(nadv_counter+1))
    print("update bld/configure: $chem_nadv = " + str(nadv_counter+1))
    print("update chem_mods.F90: $nslvd = " + str(nonadv_counter-1) + " (non-advect minus CO2)")
    print("update .xml compset files with correct dep_data_file for dry/wetdep list updates")
    print("note aerosol dry and wetdep are currently in gas list for backward compatibility")

elif mode == "wrfgc" or mode == "wrf":

    # WRF format... lots to generate.
    # registry.chem
    registry_chem_lines = "state   real   -           ikjftb   chem        1         -   -                             -\n"
    registry_chem_gas = ""
    registry_chem_nongas = ""

    # gigc_set_wrf (GC -> WRF, update all entries)
    gigc_set_wrf_other = ""
    gigc_set_wrf_simplesoa = ""
    gigc_set_wrf_complexsoa = ""
    gigc_set_wrf_complexsoa_svpoa = ""

    # gigc_get_wrf (WRF -> GC, only advected)
    gigc_get_wrf_other = ""
    gigc_get_wrf_simplesoa = ""
    gigc_get_wrf_complexsoa = ""
    gigc_get_wrf_complexsoa_svpoa = ""

    # idxarray
    gigc_idxarray = "integer :: "

    # lookup code
    gigc_idxlookup = ""

    # debug output
    gigc_idxdebug = ""

    # mozbc .inp file
    mozbc_input_file = """&control

do_bc     = .true.,
do_ic     = .true.,
domain    = 1,

dir_wrf = '/home/hplin/wrfgc/WRF/run/'
dir_moz = '/home/hplin/wrfgc/mozbc/'
fn_moz  = 'wrfgc_icbc_data_from_matlab.nc'\n\nspc_map ="""

    # for WRF-GC, for prettiness, the parsed_spc_list is sorted first
    parsed_spc_list.sort(key = lambda x: x[idx_spcName])

    for spc_idx, spc in enumerate(parsed_spc_list):
        wrfName = spc[idx_spcName].lower()
        gcName = spc[idx_spcName]

        spcName_quoted = "\"" + wrfName + "\""
        spcName_desc_quoted = "\"" + gcName + " concentration\""

        if is_gc_14_and_above:
            # use new State_Chm%Species(id_Hg2)%Conc(:,:,:) format
            get_spec = "State_Chm%Species(gi_" + wrfName + ")%Conc(II, JJ, k) = chem(i, k, j, p_" + wrfName + ") * 1.0e-6_fp\n"
            set_spec = "chem(i, k, j, p_" + wrfName + ") = State_Chm%Species(gi_" + wrfName + ")%Conc(II, JJ, k) * 1.0e+6_fp\n"
        else:
            get_spec = "State_Chm%Species(II, JJ, k, gi_" + wrfName + ") = chem(i, k, j, p_" + wrfName + ") * 1.0e-6_fp\n"
            set_spec = "chem(i, k, j, p_" + wrfName + ") = State_Chm%Species(II, JJ, k, gi_" + wrfName + ") * 1.0e+6_fp\n"

        if spc[idx_Advect] == "T":
            registry_chem_lines += "state   real   " + wrfName.ljust(10) + "  ikjftb   chem        1         -   i0{12}rhusdf=(bdy_interp:dt)  " + spcName_quoted.ljust(14) + " " + spcName_desc_quoted.ljust(27) + " \"ppmv\"\n"

            # add to gigc_get_wrf
            if gcName in complexsoa:
                gigc_get_wrf_complexsoa += "    " + get_spec
            elif gcName in complexsoa_svpoa:
                gigc_get_wrf_complexsoa_svpoa += "    " + get_spec
            elif gcName in simplesoa:
                gigc_get_wrf_simplesoa += "    " + get_spec
            else:
                gigc_get_wrf_other += get_spec

            # add to mozbc spec
            mozbc_input_file += "        '" + wrfName + " -> " + gcName + "',\n"

        else:
            registry_chem_lines += "state   real   " + wrfName.ljust(10) + "  ikjft    chem        1         -         rhusdf=(bdy_interp:dt)  " + spcName_quoted.ljust(14) + " " + spcName_desc_quoted.ljust(27) + " \"ppmv\"\n"

        # check if gas / nongas
        if spc[idx_spcName] in nongas:
            # non-gas
            registry_chem_nongas += wrfName + ","
        else:
            # gas
            registry_chem_gas += wrfName + ","

        # add to idx lookup
        gigc_idxarray    += "gi_" + wrfName + ","
        gigc_idxlookup   += "gi_" + wrfName + " = IND_('" + gcName + "')\n"
        gigc_idxdebug    += "write(6,*) p_" + wrfName + ", \"" + wrfName + " = " + gcName + "\", gi_" + wrfName + "\n"


        # add to gigc_set_wrf
        if gcName in complexsoa:
            gigc_set_wrf_complexsoa += "    " + set_spec
        elif gcName in complexsoa_svpoa:
            gigc_get_wrf_complexsoa_svpoa += "    " + set_spec
        elif gcName in simplesoa:
            gigc_set_wrf_simplesoa += "    " + set_spec
        else:
            gigc_set_wrf_other += set_spec

    # Assemble the mozbc file
    mozbc_input_file += "\n/\n"

    # Do some post-processing of the complexSOA code ...
    gigc_set_wrf_complexsoa       = "if(Input_Opt%LSOA) then\n" + gigc_set_wrf_complexsoa + "else\n" + gigc_set_wrf_simplesoa + "\nendif\n"
    gigc_set_wrf_complexsoa_svpoa = "if(Input_Opt%LSVPOA) then\n" + gigc_set_wrf_complexsoa_svpoa + "endif\n"

    gigc_get_wrf_complexsoa       = "if(Input_Opt%LSOA) then\n" + gigc_get_wrf_complexsoa + "else\n" + gigc_get_wrf_simplesoa + "\nendif\n"
    gigc_get_wrf_complexsoa_svpoa = "if(Input_Opt%LSVPOA) then\n" + gigc_get_wrf_complexsoa_svpoa + "endif\n"

    # Assemble the final registry entry
    registry_chem_nongas = "".join(registry_chem_nongas.rsplit(",", 1))
    registry_chem_gas    = "".join(registry_chem_gas.rsplit(",", 1))
    gigc_idxarray        = "".join(gigc_idxarray.rsplit(",", 1))

    registry_chem = "package   gchp                  chem_opt==233       -             chem:" + registry_chem_gas + "," + registry_chem_nongas + "," + wrf_extra_coupling

    # Assemble the final get routines in loop...
    gigc_get_wrf = gigc_get_wrf_other + "\n\n! Complex SOA species, only if available (hplin, 8/11/22)\n" + gigc_get_wrf_complexsoa + "\n" + gigc_get_wrf_complexsoa_svpoa
    gigc_set_wrf = gigc_set_wrf_other + "\n\n! Complex SOA species, only if available (hplin, 8/11/22)\n" + gigc_set_wrf_complexsoa + "\n" + gigc_set_wrf_complexsoa_svpoa

    # write all these to files...
    f = open("out/wrfgc_convert_state_mod_get.F90", "w")
    f.write(gigc_get_wrf)
    f.close()

    f = open("out/wrfgc_convert_state_mod_set.F90", "w")
    f.write(gigc_set_wrf)
    f.close()

    f = open("out/wrfgc_convert_state_mod_idxsetup.F90", "w")
    f.write(gigc_idxlookup + "\n\n" + gigc_idxdebug)
    f.close()

    f = open("out/wrfgc_convert_state_mod_idxdef.F90", "w")
    f.write(gigc_idxarray)
    f.close()

    f = open("out/registry.chem_wrfgc_entries.txt", "w")
    f.write(registry_chem_lines)
    f.close()

    f = open("out/mozbc.inp", "w")
    f.write(mozbc_input_file)
    f.close()

    print("---- registry.chem chem_opt = 233 entry: ----")
    print(registry_chem)
    print("---------------------------------------------")
    print("Update get_last_gas in module_input_chem_data.F to correspond with the above registry entry.\n")
    print("Update get_last_gas in module_chem_share.cpy (WRFv4+) to correspond with the above registry entry.\n")

else:
    raise ValueError("Unrecognized mode option! Set wrfgc or cesm.")
