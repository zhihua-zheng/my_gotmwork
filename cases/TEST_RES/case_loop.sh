# This is the main loop over turbulent methods, vertical resolution
# and time step.
# Used by: OCSPapa, OSMOSIS
#
# Qing Li, 20180507

#######################################################################
#                           Loop over cases                           #
#######################################################################

# name of the turbulence model
turblist=(KPP-CVMix KPPLT-EFACTOR KPPLT-ENTR OSMOSIS EPBL SMC SMCLT)

# vertical resolution
#  1 m
#  5 m
#  typical regional models (e.g., ROMS)
#  typical GCMs (e.g., CESM)
vrlist=(1m 5m)

# time step
#  1 min
#  30 min
dtlist=(60 1800)

# output file name
outname="gotm_out"

# loop over turbulent methods
for turbmethod in ${turblist[@]}; do

# loop over vertical resolution
for vr in ${vrlist[@]}; do
    case ${vr} in
        "1m")
            grid_method=0
            ddu=0
            ddl=0
            let nlev=depth
            ;;
        "5m")
            grid_method=0
            ddu=0
            ddl=0
            let nlev=depth/5
            ;;
        *)
            echo "Vertical resolution ${vr} not supported. Stop."
            exit 1
    esac

# loop over time step
for dt in ${dtlist[@]}; do

    # case name
    casename="TEST_RES/${title}/${turbmethod}_VR${vr}_DT${dt}s"
    echo ${casename}

    # set output frequency (3-hourly output)
    let nsave=10800/dt

    # create run directory
    rundir="${scratchdir}/${casename}"
    mkdir -p ${rundir}
    cd ${rundir}

    # copy base case
    cp ${basecase}/*.dat ./
    cp ${basecase}/*.nml ./

    if [ -f ${basecase}/*.gz ]; then
        cp ${basecase}/*.gz ./
        # decompress input data
        for f in *.gz; do
            gunzip -f ${f}
        done
    fi

    # set run parameters
    if [ -n "${datestart}" ] && [ -n "${dateend}" ]; then
        start_time="${datestart:0:4}-${datestart:4:2}-${datestart:6:2} 00:00:00"
        stop_time="${dateend:0:4}-${dateend:4:2}-${dateend:6:2} 00:00:00"
        ${cmd_nmlchange} -f gotmrun.nml -e start -v "${start_time}"
        ${cmd_nmlchange} -f gotmrun.nml -e stop -v "${stop_time}"
    fi
    ${cmd_nmlchange} -f gotmrun.nml -e title -v ${title}
    ${cmd_nmlchange} -f gotmrun.nml -e out_fn -v ${outname}
    ${cmd_nmlchange} -f gotmrun.nml -e dt -v ${dt}
    ${cmd_nmlchange} -f gotmrun.nml -e nsave -v ${nsave}
    ${cmd_nmlchange} -f gotmrun.nml -e nlev -v ${nlev}
    ${cmd_nmlchange} -f gotmmean.nml -e grid_method -v ${grid_method}
    ${cmd_nmlchange} -f gotmmean.nml -e ddu -v ${ddu}
    ${cmd_nmlchange} -f gotmmean.nml -e ddl -v ${ddl}

    # set turbulence method
    source ${scpt_case_turbmethod}

    # run
    gotm 2> log.${outname}

    # clean up input data
    if [ $? == 0 ]; then
        rm -f *.dat
    fi

done
done
done