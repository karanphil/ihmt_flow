#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    bindings = [
        "ihmt_num_threads":"$params.ihmt_num_threads",
        "bet_T1w":"$params.bet_T1w",
        "single_echo":"$params.single_echo",
        "filtering":"$params.filtering"]

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
    .fromPath("$params.input/**/*ihmt.nii.gz", maxDepth: 2)
    .map { [it.parent.name, it] }
    .groupTuple()
    .set{ ihmt_files_for_coregistration }

Channel
    .fromPath("$params.input/**/*.json", maxDepth: 2)
    .map { [it.parent.name, it] }
    .groupTuple()
    .set{ jsons_for_ihmt }

ref_t1 = Channel.fromFilePairs(
        "$params.input/**/*{ref.nii.gz,t1_brain_on_b0.nii.gz}",
        maxDepth: 1, size: 2, flat: true) { it.parent.name }

ihmt_files_for_coregistration
    .join(jsons_for_ihmt)
    .join(ref_t1)
    .into{ ihmt_ref_for_coregister;ihmt_ref_tmp_for_coregister }

Channel
    .fromPath("$params.input/**/*b1.nii.gz", maxDepth: 2)
    .map { [it.parent.name, it] }
    .groupTuple()
    .into{ b1_for_ihmt;check_b1 }

ihmt_ref_tmp_for_coregister
    .join(b1_for_ihmt)
    .set{ ihmt_ref_b1_for_coregister }

check_b1.count().set{ b1_counter }

process Compute_ihMT_B1 {
    cpus params.ihmt_num_threads

    input:
    set sid, file(ihmt_images), file(ihmt_json), file(ref),
        file(t1_on_b0), file(b1) from ihmt_ref_b1_for_coregister
    val(b1_count) from b1_counter

    output:
    file "${sid}__ihMTsat_b1_corrected_warped.nii.gz"
    file "${sid}__ihMTR_b1_corrected_warped.nii.gz"
    file "${sid}__MTsat_b1_corrected_warped.nii.gz"
    file "${sid}__MTR_b1_corrected_warped.nii.gz"

    when:
    b1_count > 0

    shell:
    single_echo=params.single_echo ? '--single_echo ' : ''
    filtering=params.filtering ? '--filtering ' : ''

    '''
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=!params.ihmt_num_threads
    mkdir coregistered_images
    mkdir bet_images
    mkdir segmentation
    mv !{ihmt_json} bet_images

    for image in !{ihmt_images}
    do
        basename=$(basename $image .nii.gz)
        antsRegistrationSyN.sh -d 3 -n !params.ihmt_num_threads -t r -f !{ref} \
            -m $image -o coregistered_images/${basename}
    done

    base_name_b1=$(basename !{b1} .nii.gz)
    antsApplyTransforms -d 3 -i !{b1} -r !{ref}\
        -t coregistered_images/*echo-1_acq-T1w_ihmt0GenericAffine.mat\
        -o coregistered_images/${base_name_b1}Warped.nii.gz -v

    bet coregistered_images/*echo-1_acq-T1w_ihmtWarped.nii.gz\
        coregistered_images/!{sid}__T1w_ihmtWarped_bet.nii.gz -m -R -f !{params.bet_T1w}
    mv coregistered_images/!{sid}__T1w_ihmtWarped_bet_mask.nii.gz\
        bet_images/!{sid}__brain_mask.nii.gz

    for image in coregistered_images/*ihmtWarped.nii.gz
    do
        base_name=$(basename $image Warped.nii.gz)
        ImageMath 3 bet_images/${base_name}.nii.gz\
            m $image bet_images/*brain_mask.nii.gz
    done

    base_name_b1=$(basename !{b1} Warped.nii.gz)
    ImageMath 3 bet_images/${base_name_b1}.nii.gz\
        m coregistered_images/*b1Warped.nii.gz bet_images/*brain_mask.nii.gz

    antsAtroposN4.sh -d 3\
                     -u 0\
                     -a bet_images/*echo-1_acq-T1w_ihmt.nii.gz\
                     -x bet_images/*brain_mask.nii.gz\
                     -c 3\
                     -o segmentation/!{sid}

    ImageMath 3 segmentation/tmp.nii.gz +\
        segmentation/!{sid}SegmentationPosteriors1.nii.gz\
        segmentation/!{sid}SegmentationPosteriors2.nii.gz

    ImageMath 3 segmentation/!{sid}__concatenated_mask.nii.gz +\
        segmentation/!{sid}SegmentationPosteriors3.nii.gz segmentation/tmp.nii.gz

    mrcalc segmentation/!{sid}__concatenated_mask.nii.gz 0.5 -ge\
        segmentation/!{sid}__concatenated_mask.nii.gz -force

    ImageMath 3 segmentation/!{sid}__concatenated_mask.nii.gz GE\
        segmentation/!{sid}__concatenated_mask.nii.gz 1

    scil_image_math.py convert segmentation/!{sid}__concatenated_mask.nii.gz\
        segmentation/!{sid}__concatenated_mask.nii.gz --data_type int8 -f

    scil_compute_ihMT_maps.py . segmentation/!{sid}__concatenated_mask.nii.gz\
        --in_altnp bet_images/*altnp*.nii.gz --in_altpn bet_images/*altpn*.nii.gz\
        --in_mtoff bet_images/*mtoff*.nii.gz --in_negative bet_images/*neg*.nii.gz\
        --in_positive bet_images/*pos*.nii.gz --in_t1w bet_images/*T1w*.nii.gz\
        --in_B1_map bet_images/*b1*.nii.gz --out_prefix !{sid}_ !{single_echo} !{filtering}

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
    mv outputWarped.nii.gz !{sid}__MTsat_b1_corrected_warped.nii.gz
    mv output0GenericAffine.mat !{sid}__output0GenericAffine.mat
    mv output1InverseWarp.nii.gz !{sid}__output1InverseWarp.nii.gz
    mv output1Warp.nii.gz !{sid}__output1Warp.nii.gz

    antsApplyTransforms -d 3 -i ihMT_native_maps/*_MTR*.nii.gz\
        -r !{sid}__MTsat_b1_corrected_warped.nii.gz \
        -o !{sid}__MTR_b1_corrected_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i ihMT_native_maps/*ihMTR*.nii.gz\
        -r !{sid}__MTsat_b1_corrected_warped.nii.gz \
        -o !{sid}__ihMTR_b1_corrected_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

     antsApplyTransforms -d 3 -i ihMT_native_maps/*ihMTsat*.nii.gz\
        -r !{sid}__MTsat_b1_corrected_warped.nii.gz \
        -o !{sid}__ihMTsat_b1_corrected_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat
    '''
}

process Compute_ihMT {
    cpus params.ihmt_num_threads

    input:
    set sid, file(ihmt_images), file(ihmt_json), file(ref),
        file(t1_on_b0) from ihmt_ref_for_coregister
    val(b1_count) from b1_counter

    output:
    file "${sid}__ihMTsat_warped.nii.gz"
    file "${sid}__ihMTR_warped.nii.gz"
    file "${sid}__MTsat_warped.nii.gz"
    file "${sid}__MTR_warped.nii.gz"

    when:
    b1_count == 0

    shell:
    single_echo=params.single_echo ? '--single_echo ' : ''
    filtering=params.filtering ? '--filtering ' : ''

    '''
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=!params.ihmt_num_threads
    mkdir coregistered_images
    mkdir bet_images
    mkdir segmentation
    mv !{ihmt_json} bet_images

    for image in !{ihmt_images}
    do
        basename=$(basename $image .nii.gz)
        antsRegistrationSyN.sh -d 3 -n !params.ihmt_num_threads -t r -f !{ref} \
            -m $image -o coregistered_images/${basename}
    done

    bet coregistered_images/*echo-1_acq-T1w_ihmtWarped.nii.gz\
        coregistered_images/!{sid}__T1w_ihmtWarped_bet.nii.gz -m -R -f !{params.bet_T1w}
    mv coregistered_images/!{sid}__T1w_ihmtWarped_bet_mask.nii.gz\
        bet_images/!{sid}__brain_mask.nii.gz

    for image in coregistered_images/*ihmtWarped.nii.gz
    do
        base_name=$(basename $image Warped.nii.gz)
        ImageMath 3 bet_images/${base_name}.nii.gz\
            m $image bet_images/*brain_mask.nii.gz
    done

    antsAtroposN4.sh -d 3\
                     -u 0\
                     -a bet_images/*echo-1_acq-T1w_ihmt.nii.gz\
                     -x bet_images/*brain_mask.nii.gz\
                     -c 3\
                     -o segmentation/!{sid}

    ImageMath 3 segmentation/tmp.nii.gz +\
        segmentation/!{sid}SegmentationPosteriors1.nii.gz\
        segmentation/!{sid}SegmentationPosteriors2.nii.gz

    ImageMath 3 segmentation/!{sid}__concatenated_mask.nii.gz +\
        segmentation/!{sid}SegmentationPosteriors3.nii.gz segmentation/tmp.nii.gz

    mrcalc segmentation/!{sid}__concatenated_mask.nii.gz 0.5 -ge\
        segmentation/!{sid}__concatenated_mask.nii.gz -force

    ImageMath 3 segmentation/!{sid}__concatenated_mask.nii.gz GE\
        segmentation/!{sid}__concatenated_mask.nii.gz 1

    scil_image_math.py convert segmentation/!{sid}__concatenated_mask.nii.gz\
        segmentation/!{sid}__concatenated_mask.nii.gz --data_type int8 -f

    scil_compute_ihMT_maps.py . segmentation/!{sid}__concatenated_mask.nii.gz\
        --in_altnp bet_images/*altnp*.nii.gz --in_altpn bet_images/*altpn*.nii.gz\
        --in_mtoff bet_images/*mtoff*.nii.gz --in_negative bet_images/*neg*.nii.gz\
        --in_positive bet_images/*pos*.nii.gz --in_t1w bet_images/*T1w*.nii.gz\
        --out_prefix !{sid}_ !{single_echo} !{filtering}

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
    mv outputWarped.nii.gz !{sid}__MTsat_warped.nii.gz
    mv output0GenericAffine.mat !{sid}__output0GenericAffine.mat
    mv output1InverseWarp.nii.gz !{sid}__output1InverseWarp.nii.gz
    mv output1Warp.nii.gz !{sid}__output1Warp.nii.gz

    antsApplyTransforms -d 3 -i ihMT_native_maps/*_MTR*.nii.gz\
        -r !{sid}__MTsat_warped.nii.gz \
        -o !{sid}__MTR_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

    antsApplyTransforms -d 3 -i ihMT_native_maps/*ihMTR*.nii.gz\
        -r !{sid}__MTsat_warped.nii.gz \
        -o !{sid}__ihMTR_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat

     antsApplyTransforms -d 3 -i ihMT_native_maps/*ihMTsat*.nii.gz\
        -r !{sid}__MTsat_warped.nii.gz \
        -o !{sid}__ihMTsat_warped.nii.gz -n Linear \
        -t !{sid}__output1Warp.nii.gz !{sid}__output0GenericAffine.mat
    '''
}
