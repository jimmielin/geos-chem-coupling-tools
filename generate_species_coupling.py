############################################
# Generate species coupling code.
# Used for GEOS-Chem with CESM and WRF-GC
# (c) 2022 Haipeng Lin <myself@jimmielin.me>
#
# Licensed under the GNU General Public License v2
############################################

import re
from operator import itemgetter

# choose parser_mode: "regex" for 13.4.x, or "yaml" for 14.1.0+
parser_mode = "yaml"

# choose input file:
input_file_name = "v14.0.0.yml"

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

                # heap alloc.
                parsed_spc_list.append((N, spcName, Is_Advected, Is_DryDep, Is_WetDep, Mw_g))

            print(parsed_spc_list)
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

# cesm | wrfgc
mode = "cesm"

# sort by adv to keep them first - even though GC should keep them first
# this is for sanity.
parsed_spc_list.sort(key = lambda x: 0 if x[idx_Advect] == "T" else 1)

#
# These are "quirky" species, which necessitate special handling
# in coupled models such as WRF and CESM.
#

# Non-gas-phase species ("Aerosols")
nongas = ["AERI", "ASOAN", "ASOA1", "ASOA2", "ASOA3", "ASOG1", "ASOG2", "ASOG3", "AONITA", "BCPI", "BCPO", "BrSALA", "BrSALC", "DMS", "DST1", "DST2", "DST3", "DST4", "INDIOL", "IONITA", "ISALA", "ISALC", "ISN1OA", "ISN1OG", "ISOA1", "ISOA2", "ISOA3", "MONITA", "MSA", "NH4", "NIT", "NITs", "OCPI", "OCPO", "OPOA1", "OPOA2", "POA1", "POA2", "SALA", "SALC", "SALACl", "SALAAL", "SALCAL", "SALCCl", "SO4", "SO4s", "SOAIE", "SOAGX", "SOAME", "SOAMG", "SOAP", "SOAS", "TSOA0", "TSOA1", "TSOA2", "TSOA3", "TSOG0", "TSOG1", "TSOG2", "TSOG3", "pFe", "POA1", "POA2", "POG1", "POG2"]

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

# WRF also requires mapping of species to WRF-Chem names for
# compatibility. TODO hplin 7/20/22


############ NO USER CONFIGURABLE CODE BELOW ############
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
    first_nonadv = True

    # Final strings. Headers will be filled later, as we need to count.
    solsym = ""
    adv_mass = ""
    drydep_list = ""
    wetdep_list = ""

    for spc_idx, spc in enumerate(parsed_spc_list):
        # invariants do not need to be skipped, they are actually duplicated
        # if spc[1] in cesm_invariants:
        #     continue

        # insert drydep / wetdep
        if spc[idx_DryDep] == "T" and (not spc[idx_spcName].upper() in nongas_upper):
            drydep_list += "'" + spc[idx_spcName].upper() + "',"
        if spc[idx_WetDep] == "T" and (not spc[idx_spcName].upper() in nongas_upper):
            wetdep_list += "'" + spc[idx_spcName].upper() + "',"
        
        # insert aerosols in between T and F. check for first F
        if first_nonadv and spc[2] == "F":
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
                    solsym += "&\r\n"
                    adv_mass += "&\r\n"
        
        solsym += "'" + spc[idx_spcName].ljust(15) + "', "
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

    print(solsym)
    print(adv_mass)
    print("<drydep_list>" + drydep_list + "</drydep_list>")
    print("<gas_wetdep_list>" + wetdep_list + "</gas_wetdep_list>")
    print("update chem_mods.F90: gas_pcnst = " + str(total_counter))
    print("update chem_mods.F90: nTracersMax = " + str(nadv_counter+1))
    print("update bld/configure: $chem_nadv = " + str(nadv_counter+1))