#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    bindings = [
        "ihmt_num_threads":"$params.ihmt_num_threads",
        "bet_T1w":"$params.bet_T1w",
        "filtering":"$params.filtering",
        "extended":"$params.extended"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "SCIL ihMT pipeline"
log.info "=================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "Execution time: $workflow.duration"
}

Channel
    .fromPath("$params.input/**/*ihmt.nii.gz", maxDepth: 1)
    .map { [it.parent.name, it] }
    .groupTuple()
    .into{ ihmt_files_for_coregistration; check_single_multi_echo; check_mix }

Channel
    .fromPath("$params.input/**/*.json", maxDepth: 1)
    .map { [it.parent.name, it] }
    .groupTuple()
    .set{ jsons_for_ihmt }

ref_t1 = Channel.fromFilePairs(
        "$params.input/**/*{ref.nii.gz,t1_brain_on_b0.nii.gz}",
        maxDepth: 1, size: 2, flat: true) { it.parent.name }

ihmt_files_for_coregistration
    .join(jsons_for_ihmt)
    .join(ref_t1)
    .into{ ihmt_ref_for_coregister;ihmt_ref_for_count }

Channel
    .fromPath("$params.input/**/*b1.nii.gz", maxDepth: 1)
    .map { [it.parent.name, it] }
    .groupTuple()
    .into{ b1_for_ihmt;check_b1}

ihmt_ref_for_coregister
    .join(b1_for_ihmt, remainder:true)
    .set{ ihmt_ref_b1_for_coregister }

check_b1.count().set{ b1_counter }
ihmt_ref_for_count.count().into{number_subj_for_compare; number_of_subj_for_echo}

Channel
    .fromPath("$params.input/**/fitValues*.mat", maxDepth: 1)
    .map { [it.parent.name, it] }
    .groupTuple()
    .set{ b1_fitvalues }

ihmt_ref_b1_for_coregister
    .join(b1_fitvalues)
    .set{ ihmt_ref_b1_for_coregister }

number_subj_for_compare
    .concat(b1_counter)
    .toList()
    .subscribe{a, b -> if (a != b && b > 0)
    error "Error ~ Some subjects have a b1 image and others don't.\n" +
          "Please be sure to have the same acquisitions for all subjects."}

multi_echo = true
check_single_multi_echo.map { it[1] }
    .flatten()
    .count()
    .concat(number_of_subj_for_echo)
    .toList()
    .subscribe{a, b -> if (a % 6 == 0 && b * 6 == a)
    multi_echo = false}

check_mix.map { it[1].size() }
    .unique()
    .toList()
    .subscribe{ List a -> if (a.size() > 1)
    error "Error ~ Some subjects have a multi-echo acquisitions and others don't.\n" +
          "Please be sure to have the same acquisitions for all subjects."}

number_subj_for_compare
    .concat(b1_counter)
    .toList()
    .subscribe{a, b -> if (a == b)
    println "-> Use B1 images for correction\n"}

process Compute_ihMT {
    cpus params.ihmt_num_threads

    input:
    set sid, file(ihmt_images), file(ihmt_json), file(ref),
        file(t1_on_b0), file(b1), file(fitvalues) from ihmt_ref_b1_for_coregister
    val(b1_count) from b1_counter

    output:
    file("Register_ihMT_maps")
    file("Register_contrast_maps")
    file("Register_MT_T1w")

    file("Contrasts_ihMT_maps") optional true
    file("ihMT_native_maps") optional true
    file("Segmentation") optional true
    file("Coregistered_images") optional true
    file("Bet_images") optional true
    file("Registration_files") optional true

    shell:
    filtering=params.filtering ? '--filtering ' : ''
    b1_input_params=b1_count ? '--in_B1_map Bet_images/*b1*.nii.gz' : ''
    b1_method_params=b1_count ? '--B1_correction_method model_based' : ''
    b1_fitvalues_params=b1_count ? '--in_B1_fitvalues B1_fitValues/fitValues_SP_1.mat B1_fitValues/fitValues_SN_1.mat B1_fitValues/fitValues_D_1.mat ' : ''
    b1_ext=b1_count ? '_b1' : ''
    extended=params.extended ? 'true' : 'false'

    '''    
    echo !{multi_echo}
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=!{params.ihmt_num_threads}
    
    mkdir Contrasts_ihMT_maps
    mkdir ihMT_native_maps
    mkdir Segmentation
    mkdir Coregistered_images
    mkdir Bet_images
    mkdir Registration_files

    mkdir Register_ihMT_maps
    mkdir Register_contrast_maps
    mkdir Register_MT_T1w
    
    mv !{ihmt_json} Bet_images

    mkdir B1_fitValues
    mv !{fitvalues} B1_fitValues

    for image in !{ihmt_images}
    do
        basename=$(basename $image .nii.gz)
        antsRegistrationSyN.sh -d 3 -n !{params.ihmt_num_threads} -t r -f !{ref} \
            -m $image -o Coregistered_images/${basename}
    done

    if [[ !{b1_count} != 0 ]]
    then
        base_name_b1=$(basename !{b1} .nii.gz)
        antsApplyTransforms -d 3 -i !{b1} -r !{ref}\
            -t Coregistered_images/*echo-1_acq-T1w_ihmt0GenericAffine.mat\
            -o Coregistered_images/${base_name_b1}Warped.nii.gz -n Linear -v
    fi

    bet Coregistered_images/*echo-1_acq-T1w_ihmtWarped.nii.gz\
        Coregistered_images/!{sid}__T1w_ihmtWarped_bet.nii.gz -m -R -f !{params.bet_T1w}
    mv Coregistered_images/!{sid}__T1w_ihmtWarped_bet_mask.nii.gz\
        Bet_images/!{sid}__brain_mask.nii.gz

    for image in Coregistered_images/*ihmtWarped.nii.gz
    do
        base_name=$(basename $image Warped.nii.gz)
        ImageMath 3 Bet_images/${base_name}.nii.gz\
            m $image Bet_images/*brain_mask.nii.gz
    done

    if [[ !{b1_count} != 0 ]]
    then
        base_name_b1=$(basename !{b1} Warped.nii.gz)
        ImageMath 3 Bet_images/${base_name_b1}.nii.gz\
            m Coregistered_images/*b1Warped.nii.gz Bet_images/*brain_mask.nii.gz
    fi

    antsAtroposN4.sh -d 3\
                     -u 0\
                     -a Bet_images/*echo-1_acq-T1w_ihmt.nii.gz\
                     -x Bet_images/*brain_mask.nii.gz\
                     -c 3\
                     -o Segmentation/!{sid}

    ImageMath 3 Segmentation/tmp.nii.gz +\
        Segmentation/!{sid}SegmentationPosteriors1.nii.gz\
        Segmentation/!{sid}SegmentationPosteriors2.nii.gz

    ImageMath 3 Segmentation/!{sid}__concatenated_mask.nii.gz +\
        Segmentation/!{sid}SegmentationPosteriors3.nii.gz Segmentation/tmp.nii.gz

    mrcalc Segmentation/!{sid}__concatenated_mask.nii.gz 0.5 -ge\
        Segmentation/!{sid}__concatenated_mask.nii.gz -force

    ImageMath 3 Segmentation/!{sid}__concatenated_mask.nii.gz GE\
        Segmentation/!{sid}__concatenated_mask.nii.gz 1

    scil_image_math.py convert Segmentation/!{sid}__concatenated_mask.nii.gz\
        Segmentation/!{sid}__concatenated_mask.nii.gz --data_type int8 -f

    export single_echo="--single_echo"
    if [[ !{multi_echo} == true ]]
    then
        export single_echo=""
    fi

    scil_compute_ihMT_maps.py . Segmentation/!{sid}__concatenated_mask.nii.gz\
        --in_altnp Bet_images/*altnp*.nii.gz --in_altpn Bet_images/*altpn*.nii.gz\
        --in_mtoff Bet_images/*mtoff*.nii.gz --in_negative Bet_images/*neg*.nii.gz\
        --in_positive Bet_images/*pos*.nii.gz --in_t1w Bet_images/*T1w*.nii.gz\
        --out_prefix !{sid}_ ${single_echo} !{filtering} !{b1_input_params}\
        !{b1_method_params} !{b1_fitvalues_params}

    base_name_mt=$(basename ihMT_native_maps/*_MTsat*)

    antsRegistration --dimensionality 3 --float 0\
        --output [output,outputWarped.nii.gz,outputInverseWarped.nii.gz]\
        --interpolation Linear --use-histogram-matching 0\
        --winsorize-image-intensities [0.005,0.995]\
        --initial-moving-transform [!{t1_on_b0},ihMT_native_maps/${base_name_mt},1]\
        --transform Rigid['0.2']\
        --metric MI[!{t1_on_b0},ihMT_native_maps/${base_name_mt},1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform Affine['0.2']\
        --metric MI[!{t1_on_b0},ihMT_native_maps/${base_name_mt},1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform SyN[0.1,3,0]\
        --metric MI[!{t1_on_b0},ihMT_native_maps/${base_name_mt},1,32]\
        --convergence [50x25x10,1e-6,10] --shrink-factors 4x2x1\
        --smoothing-sigmas 3x2x1
    mv outputWarped.nii.gz !{sid}__MTsat!{b1_ext}_warped.nii.gz
    mv output0GenericAffine.mat !{sid}__output0GenericAffine.mat
    mv output1InverseWarp.nii.gz !{sid}__output1InverseWarp.nii.gz
    mv output1Warp.nii.gz !{sid}__output1Warp.nii.gz

    antsApplyTransforms -d 3 -i ihMT_native_maps/*_MTR*.nii.gz\
        -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
        -o Register_ihMT_maps/!{sid}__MTR!{b1_ext}_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i ihMT_native_maps/*ihMTR*.nii.gz\
        -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
        -o Register_ihMT_maps/!{sid}__ihMTR!{b1_ext}_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

     antsApplyTransforms -d 3 -i ihMT_native_maps/*ihMTsat*.nii.gz\
        -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
        -o Register_ihMT_maps/!{sid}__ihMTsat!{b1_ext}_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i Contrasts_ihMT_maps/*_positive*.nii.gz\
        -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
        -o Register_contrast_maps/!{sid}__positive!{b1_ext}_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i Contrasts_ihMT_maps/*negative*.nii.gz\
        -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
        -o Register_contrast_maps/!{sid}__negative!{b1_ext}_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i Contrasts_ihMT_maps/*altnp*.nii.gz\
        -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
        -o Register_contrast_maps/!{sid}__altnp!{b1_ext}_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i Contrasts_ihMT_maps/*altpn*.nii.gz\
       -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
       -o Register_contrast_maps/!{sid}__altpn!{b1_ext}_warped.nii.gz -n Linear \
       -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i Contrasts_ihMT_maps/*reference*.nii.gz\
       -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
       -o Register_contrast_maps/!{sid}__reference!{b1_ext}_warped.nii.gz -n Linear \
       -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

   antsApplyTransforms -d 3 -i Contrasts_ihMT_maps/*T1w*.nii.gz\
      -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
      -o Register_MT_T1w/!{sid}__T1w!{b1_ext}_warped.nii.gz -n Linear \
      -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    if [[ !{b1_count} != 0 ]]
    then
        antsApplyTransforms -d 3 -i Bet_images/*b1*.nii.gz\
            -r !{sid}__MTsat!{b1_ext}_warped.nii.gz \
            -o Register_contrast_maps/!{sid}_!{b1_ext}_warped.nii.gz -n Linear \
            -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat
    fi


    mv !{sid}__MTsat!{b1_ext}_warped.nii.gz Register_ihMT_maps/
    mv !{sid}__output1Warp.nii.gz Registration_files/
    mv !{sid}__output0GenericAffine.mat Registration_files/

    if [[ !{extended} == "false" ]]
    then
    mv Contrasts_ihMT_maps Contrasts_ihMT_maps_tmp 
    mv ihMT_native_maps ihMT_native_maps_tmp
    mv Segmentation Segmentation_tmp
    mv Coregistered_images Coregistered_images_tmp
    mv Bet_images Bet_images_tmp
    mv Registration_files Registration_files_tmp
    fi
    '''
}
