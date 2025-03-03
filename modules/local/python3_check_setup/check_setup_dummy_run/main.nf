#!/usr/bin/env nextflow

// Don't actually perfomr any checks, but this process is required for later steps to complete
process CHECK_SETUP_DUMMY_RUN {

    output:
        val('CONTROL_1')    // Prevents later processes executing until after the check is completed

    script:
    """
    echo No checks needed
    """
}