# AutoFitLor.tcl
#
# This script automates the fitting of Lorentzian models to spectral data in XSPEC.
# It iteratively adds Lorentzian components to fit power density spectra (PDS),
# real and imaginary parts, phase lags, and coherence spectra, using statistical tests
# to determine the optimal number of Lorentzians.
#

# Add the Tcllib path for math::statistics package
# CHANGE THIS TO OWN PATH or EDIT PATH BEHOREHAND AND COMMENT 
lappend auto_path /opt/software/tcllib/
package require math::statistics

# -------------------
# Global Variables
# -------------------
set ::pi 3.14159265
set ::pshift [expr $::pi / 4.0]

# -------------------
# Utility Functions
# -------------------

# sigma_to_pvalue --
#   Converts a sigma value to a two-tailed p-value using the normal distribution.
# Arguments:
#   sigma - The sigma value (e.g., 3.0 for 3-sigma)
# Returns:
#   p-value as a double
proc sigma2pvalue {sigma} {
    set p [expr {2 * (1 - [::math::statistics::cdf-normal 0 1 $sigma])}]
    return $p
}


#  WW_runs_test_residuals --
#   Performs a Wald-Wolfowitz runs test on residuals to check for randomness.
#   If residuals are non-random, more Lorentzians may be needed.
# Arguments:
#   fit_data - List of spectra data, each containing [spec_name, freqs, residuals]
# Returns:
#   1 if residuals are random, 0 otherwise
proc WW_runs_test {fit_data} {
    set min_p 1.0
    set is_random 1
    foreach spec $fit_data {
        lassign $spec spec_name freqs residuals
        set binary {}
        foreach res_entry $residuals {
            lassign $res_entry res err
            lappend binary [expr {$res > 0 ? 1 : -1}]
        }
        set n_pos 0
        set n_neg 0
        set runs 1
        set prev_sign 0
        foreach val $binary {
            if {$val == 1} { incr n_pos } else { incr n_neg }
            if {$val != 0 && $prev_sign != 0 && $val != $prev_sign} { incr runs }
            if {$val != 0} { set prev_sign $val }
        }
        if {$n_pos == 0 || $n_neg == 0} {
            continue
        }
        set n [expr {$n_pos + $n_neg}]
        set expected_runs [expr {2.0 * $n_pos * $n_neg / $n + 1}]
        set variance [expr {2.0 * $n_pos * $n_neg * (2.0 * $n_pos * $n_neg - $n) / ($n * $n * ($n - 1))}]
        if {$variance == 0} {
            continue
        }
        set z [expr {($runs - $expected_runs) / sqrt($variance)}]
        set p [sigma2pvalue [expr {abs($z)}]]
        if {$p < $min_p} { set min_p $p }
        if {$p < $::pval} { set is_random 0 }
    }
    puts "WW runs test p-value: $min_p (is_random: $is_random)"
    return $is_random
}


## Durbin-Watson Test 
## ----------------------------------------
## Interpretation:
##    DW ≈ 2: No significant autocorrelation (desired for a good fit).
##    DW < 1: Positive autocorrelation (residuals show trends; suggests missing structure in the model).
##    DW > 3: Negative autocorrelation (oscillatory pattern; might indicate overfitting or incorrect error assumptions).
## How to Use:
##    Compute the Durbin-Watson statistic on your residuals.
##    If DW deviates significantly from 2, consider:
##        Adding a missing component (e.g., extra absorption or emission lines).
##        Adjusting systematics.
##        Checking for correlated noise sources.

proc DW_test {residuals} {
    set num [llength $residuals]
    if {$num < 2} {
        return "Insufficient data for Durbin-Watson test"
    }
    
    set numerator 0
    set denominator 0
    for {set i 1} {$i < $num} {incr i} {
        set diff [expr {[lindex $residuals $i] - [lindex $residuals [expr {$i-1}]]}]
        set numerator [expr {$numerator + $diff * $diff}]
    }
    
    foreach res $residuals {
        set denominator [expr {$denominator + $res * $res}]
    }
    
    return [expr {$numerator / $denominator}]
}

## AIC Test 
## -------------------------------------------------------------------
# Interpretation:
# Lower AIC values indicate a better fit (relative to other models).
# ΔAIC = AIC(model1) - AIC(model2):
#    ΔAIC > 10: Strong evidence that one model is better.
#    ΔAIC ≈ 2-6: Weak to moderate preference.
#    ΔAIC < 2: Models are nearly equivalent.

# How to Use:
#    Compute AIC for different models (e.g., simple vs. complex).
#    Compare AIC values:
#        If adding a parameter lowers AIC, it improves the fit.
#        If AIC increases, the extra parameter may not be justified.


proc AIC_test {residuals num_params} {
    set num [llength $residuals]
    if {$num == 0} {
        return "No residuals provided"
    }
    
    set sum_sq 0
    foreach res $residuals {
        set sum_sq [expr {$sum_sq + $res * $res}]
    }
    
    set log_likelihood [expr {-0.5 * $num * log($sum_sq / $num)}]
    set aic_value [expr {2 * $num_params - 2 * $log_likelihood}]
    
    return $aic_value
}


# ftest_model_improvement --
#   Performs an F-test to evaluate model improvement between pre- and post-fit.
# Arguments:
#   chi_pre  - Chi-squared value before adding a Lorentzian
#   dof_pre  - Degrees of freedom before adding a Lorentzian
#   chi_post - Chi-squared value after adding a Lorentzian
#   dof_post - Degrees of freedom after adding a Lorentzian
# Returns:
#   List of [f_stat, p_value]
proc ftest_model_improvement {chi_pre dof_pre chi_post dof_post} {
    set cmd "ftest $chi_post $dof_post $chi_pre $dof_pre"
    if {[catch {eval $cmd} err]} {
        [format "F-test failed: %s" $err]
        return [list 0.0 1.0]
    }
    set p_value [tcloutr ftest]
    set f_stat [expr {($chi_pre - $chi_post) / ($dof_pre - $dof_post) / ($chi_post / $dof_post)}]
    return [list $f_stat $p_value]
}





# get_fit_data --
#   Retrieves fitting data (frequencies, residuals, errors) for each spectrum.
# Arguments:
#   spectra - List of spectrum names (e.g., {1pds 2pds 3rea 4ima})
# Returns:
#   List of [spec_name, freqs, residuals] for each spectrum
proc get_fit_data {spectra} {

    set ::spec_map {1pds 1 2pds 2 3rea 3 4ima 4 5pla 5 6coh 6}
    set fit_data {}
    
    # First pass: Determine the maximum number of bins
    set max_bins 0
    foreach spec_name $spectra {
        set spec_num [dict get $::spec_map $spec_name]
        set freqs [lsearch -all -inline -not -exact [split [tcloutr plot data x $spec_num] " "] {}]
        set n_bins [llength $freqs]
        if {$n_bins > $max_bins} {
            set max_bins $n_bins
        }
    }

    # Second pass: Process each spectrum and pad if necessary
    foreach spec_name $spectra {
        set spec_num [dict get $::spec_map $spec_name]
        set freqs [lsearch -all -inline -not -exact [split [tcloutr plot data x $spec_num] " "] {}]
        set data [lsearch -all -inline -not -exact [split [tcloutr plot data y $spec_num] " "] {}]
        set errors [lsearch -all -inline -not -exact [split [tcloutr plot data yerr $spec_num] " "] {}]
        set model [lsearch -all -inline -not -exact [split [tcloutr plot data model $spec_num] " "] {}]
        set residuals {}
        for {set i 0} {$i < $max_bins} {incr i} {
            if {$i < [llength $data]} {
                set data_val [lindex $data $i]
                set model_val [lindex $model $i]
                set err_val [lindex $errors $i]
            } else {
                set data_val 0.0
                set model_val 0.0
                set err_val 1.0
            }
            # Validate data value
            if {$data_val eq "" || ![string is double -strict $data_val]} {
                puts "Warning: Invalid data value at index $i for spectrum $spec_name: \
                '$data_val', setting to 0"
                set data_val 0.0
            }
            # Validate model value
            if {$model_val eq "" || ![string is double -strict $model_val]} {
                puts "Warning: Invalid model value at index $i for spectrum $spec_name: \
                '$model_val', setting to 0"
                set model_val 0.0
            }
            # Compute residual
            set res [expr {$data_val - $model_val}]
            # Validate error value
            if {$err_val eq "" || ![string is double -strict $err_val] || $err_val == 0} {
                puts "Warning: Invalid error value at index $i for spectrum $spec_name: \
                '$err_val', setting to 1.0"
                set err 1.0
            } else {
                set err $err_val
            }
            lappend residuals [list $res $err]
        }
        # Pad freqs if necessary
        while {[llength $freqs] < $max_bins} {
            lappend freqs 0.0
        }
        lappend fit_data [list $spec_name $freqs $residuals]
    }
    return $fit_data
}




# find_largest_residual --
#   Identifies the largest residual peak across spectra to determine where to add a new Lorentzian.
# Arguments:
#   fit_data - List of spectra data, each containing [spec_name, freqs, residuals]
# Returns:
#   List of [peak_freq, width, norm] for the new Lorentzian
proc find_largest_residual {fit_data} {
    set freqs [lindex [lindex $fit_data 0] 1]
    set n_spectra [llength $fit_data]
    set n_points [llength $freqs]

    # Step 1: Identify runs in each spectrum
    set all_runs {}
    foreach spec $fit_data {
        lassign $spec spec_name _ residuals
        set runs {}
        set current_run {}
        set prev_sign 0
        for {set i 0} {$i < $n_points} {incr i} {
            lassign [lindex $residuals $i] res err
            set sign [expr {$res > 0 ? 1 : -1}]
            if {$i == 0 || $sign == $prev_sign} {
                lappend current_run $i
            } else {
                if {[llength $current_run] >= 3} {
                    lappend runs $current_run
                }
                set current_run [list $i]
            }
            set prev_sign $sign
        }
        if {[llength $current_run] >= 3} {
            lappend runs $current_run
        }
        lappend all_runs $runs
    }

    # Step 2: Find shared runs across spectra (First Pass)
    set shared_runs {}
    set run_id 0
    foreach run_set [lindex $all_runs 0] {
        set start_idx [lindex $run_set 0]
        set end_idx [lindex $run_set end]
        set shared 1
        set shared_spectra {0}
        for {set s 1} {$s < $n_spectra} {incr s} {
            set found 0
            foreach other_run [lindex $all_runs $s] {
                set other_start [lindex $other_run 0]
                set other_end [lindex $other_run end]
                if {$other_start <= $end_idx && $other_end >= $start_idx} {
                    set found 1
                    lappend shared_spectra $s
                    break
                }
            }
            if {!$found} {
                set shared 0
                break
            }
        }
        if {$shared} {
            lappend shared_runs [list $run_id $start_idx $end_idx $shared_spectra]
            incr run_id
        }
    }

    # Step 3: Compute total chi-squared for each shared run
    set run_scores {}
    foreach run $shared_runs {
        lassign $run run_id start_idx end_idx shared_spectra
        set total_chi2 0.0
        set max_res 0.0
        set peak_idx $start_idx
        foreach s $shared_spectra {
            set spec [lindex $fit_data $s]
            lassign $spec spec_name _ residuals
            for {set i $start_idx} {$i <= $end_idx} {incr i} {
                if {$i >= [llength $residuals]} {
                    continue
                }
                lassign [lindex $residuals $i] res err
                set chi2 [expr {($res / $err) ** 2}]
                set total_chi2 [expr {$total_chi2 + $chi2}]
                if {$res > $max_res} {
                    set max_res $res
                    set peak_idx $i
                }
            }
        }
        lappend run_scores [list $run_id $start_idx $end_idx $total_chi2 $peak_idx $max_res \
        $shared_spectra]
    }

    # Step 4: If no shared runs, look for significant runs in individual spectra (Second Pass)
    if {[llength $run_scores] == 0} {
        set run_id 0
        for {set s 0} {$s < $n_spectra} {incr s} {
            set spec [lindex $fit_data $s]
            lassign $spec spec_name _ residuals
            foreach run [lindex $all_runs $s] {
                set start_idx [lindex $run 0]
                set end_idx [lindex $run end]
                set total_chi2 0.0
                set max_res 0.0
                set peak_idx $start_idx
                for {set i $start_idx} {$i <= $end_idx} {incr i} {
                    if {$i >= [llength $residuals]} {
                        continue
                    }
                    lassign [lindex $residuals $i] res err
                    set chi2 [expr {($res / $err) ** 2}]
                    set total_chi2 [expr {$total_chi2 + $chi2}]
                    if {$res > $max_res} {
                        set max_res $res
                        set peak_idx $i
                    }
                }
                lappend run_scores [list $run_id $start_idx $end_idx $total_chi2 $peak_idx $max_res \
                [list $s]]
                incr run_id
            }
        }
    }

    # Step 5: Sort by total chi-squared (descending)
    if {[llength $run_scores] == 0} {
        return [list 1.0 0.1 1e-3]
    }
    set run_scores [lsort -real -decreasing -index 3 $run_scores]

    # Step 6: Select the most significant run
    lassign [lindex $run_scores 0] run_id start_idx end_idx total_chi2 peak_idx max_res shared_spectra
    set peak_freq [lindex $freqs $peak_idx]

    # Step 7: Compute FWHM for width
    set half_max [expr {$max_res / 2.0}]
    set left_idx $peak_idx
    set right_idx $peak_idx
    set spec [lindex $fit_data [lindex $shared_spectra 0]]
    lassign $spec spec_name _ residuals
    for {set i $peak_idx} {$i >= $start_idx} {incr i -1} {
        lassign [lindex $residuals $i] res err
        if {$res <= $half_max} {
            set left_idx $i
            break
        }
    }
    for {set i $peak_idx} {$i <= $end_idx} {incr i} {
        lassign [lindex $residuals $i] res err
        if {$res <= $half_max} {
            set right_idx $i
            break
        }
    }
    set f_left [lindex $freqs $left_idx]
    set f_right [lindex $freqs $right_idx]
    set width [expr {$f_right - $f_left}]
    if {$width <= 0} { set width [expr 0.5*$peak_freq] }

    # Step 8: Compute norm
    set norm [expr max({$max_res * $width * sqrt(2 * $::pi)},1e-4)]
    puts "Adding Lorentzian at freq=$peak_freq, width=$width, norm=$norm \
    (shared across [llength $shared_spectra] spectra, chi2=$total_chi2)"
    return [list $peak_freq $width $norm]
}


# -------------------
# Main Fitting Procedure
# -------------------

# autofitlor --
#   Main procedure to automatically fit Lorentzians to spectral data.
# Arguments:
#   basename - Base name for the output (e.g., "zone1")
#   sigma    - Sigma threshold for statistical tests (default: 5.0)
#   spectra  - List of spectra to fit (default: {1pds 2pds 3rea 4ima})
proc autofitlor {{basename "logfile"} {sigma 3.0} {spectra {1pds 2pds 3rea 4ima}}} {
    
    set ::pval [sigma2pvalue $sigma]
         
    set model_string [tcloutr model]
    if {$model_string == ""} {
    	puts "No model found. Starting plor model"
	    set ::nlor 1
	    startplor
	    fit
    } else {
    	renorm; fit
		set ::nlor [expr [tcloutr varpar]/6]
	    sort_lorentzians
    	puts "Model with ::nlor lorentzians found"
    }

   	puts "Starting fit with sigma=$sigma (pval=$::pval) threshold"
    set continue_fitting 1
    while {$continue_fitting} {

        incr ::nlor
        puts "Attempting to fit Lorentzian $::nlor"
        set fit_data [get_fit_data $spectra]
        set is_random [WW_runs_test $fit_data]
        if {$is_random} {
            incr ::nlor -1
            puts "Residuals are random, stopping at $::nlor Lorentzians"
            break
        }
        
        set ::chi_pre [tcloutr stat]
        set ::dof_pre [lindex [tcloutr dof] 0]
        lassign [find_largest_residual $fit_data] freq width norm
        
        addplor $::nlor $freq $width $norm
        renorm; fit
        shakefit
          
        set ::chi_post [tcloutr stat]
        set ::dof_post [lindex [tcloutr dof] 0]
        lassign [ftest_model_improvement $::chi_pre $::dof_pre $::chi_post $::dof_post] f_test p_value
        
        if {$p_value >= $::pval} {
            delplor $::nlor
            incr ::nlor -1
            puts "F-test rejected Lorentzian $::nlor, stopping at $::nlor Lorentzians"
            set continue_fitting 0
        }
        
    	sort_lorentzians
    }
    
    lassign [ftest_model_improvement $::chi_pre $::dof_pre $::chi_post $::dof_post] ::last_ftest ::last_pvalue 
    set ::last_ww_results [WW_runs_test [get_fit_data $spectra]]
    save_fit_log $basename $spectra $sigma
    
    
}

# -------------------
# XSPEC Model Management
# -------------------

# shake --
#   Calculates parameter errors at 5 sigma level.
proc shakefit {} {
    set m1 [tcloutr modpar 1pds]
    set m3 [tcloutr modpar 3rea]
    error 25. max 1000. 1pds:1-$m1 2pds:1-$m1 3rea:1-$m3
}


# startplor --
#   Initializes a Lorentzian model for spectra 1pds, 2pds, 3rea, 4ima, 5pla, 6coh.
proc startplor {} {
    # Check current model definitions and erase plags or cohes if present
    set model_string [tcloutr model]
    if {[string match "*xlor*" $model_string]} {
        puts "Erasing existing xlor model"
        mdefine xlor :
    }
    if {[string match "*ylor*" $model_string]} {
        puts "Erasing existing ylor model"
        mdefine ylor :
    }

    set chatlevel [scan [tcloutr chatter] "%d"]
    chatter 0
    query yes
    xset LINECRITLEVEL 1e-10
    mdefine xlor Lorentz(LineE,Width)*cos(atan2(sin(plag+$::pshift),cos(plag+$::pshift))) : add
    mdefine ylor Lorentz(LineE,Width)*sin(atan2(sin(plag+$::pshift),cos(plag+$::pshift))) : add

    puts "Starting a pLorentzian model for 1pds,2pds,3rea,4ima"
    model 1:1pds lor
    0.1  0.01 1e-10 1e-10 1e10 1e10
    0.1  0.01 1e-10 1e-10 1e10 1e10
    0.01 0.01 1e-10 1e-10 1e10 1e10

    model 2:2pds lor
    =1pds:p1
    =1pds:p2
    0.01 0.01 1e-10 1e-10 1e10 1e10

    model 3:3rea xlor
    =1pds:p1
    =1pds:p2
    0.00 0.01 -$::pi -$::pi $::pi $::pi
    0.01 0.01 1e-10 1e-10 1e10 1e10

    model 4:4ima ylor
    =3rea:1
    =3rea:2
    =3rea:3
    =3rea:4

    puts "Arranging plag and coherence to 1 pLor"
    plags 1
    cohes 1
    puts "Finished starting Lorentzian model for 1 pLor."
    chatter $chatlevel
}




# delplor --
#   Deletes a Lorentzian component from the model.
# Arguments:
#   args - Lorentzian number to delete
proc delplor {args} {
    set chatlevel [scan [tcloutr chatter] "%d"]
    chatter 0
    set nlor [lindex $args 0]
    puts "Removing Lor $nlor from 1pds,2pds,3rea,4ima"
    delco 1pds:$nlor
    delco 2pds:$nlor
    delco 3rea:$nlor
    delco 4ima:$nlor
    set lors [tcloutr modcomp 1pds]
    puts "Rearranging plag and coherence to $lors pLors"
    plags $lors
    cohes $lors
    puts "Finished deleting Lor $nlor. The model has now $lors pLors."
    chatter $chatlevel
}


# addplor --
#   Adds a Lorentzian component to the model.
# Arguments:
#   args - List of [nlor, freq, width, norm]
proc addplor {args} {
    set chatlevel [scan [tcloutr chatter] "%d"]
    chatter 0
    set nlor [lindex $args 0]
    set freq [lindex $args 1]
    set width [lindex $args 2]
    set norm [lindex $args 3]

    puts "Adding lor $nlor to 1pds"
    addco 1pds:$nlor lorentz
    $freq  0.01 1e-10 1e-10 1e10 1e10
    $width 0.01 1e-10 1e-10 1e10 1e10
    $norm  0.01 1e-10 1e-10 1e10 1e10

    puts "Adding lor $nlor to 2pds"
    addco 2pds:$nlor lorentz & /*
    newpar 2pds:[expr 3*($nlor-1)+1]=1pds:p[expr 3*($nlor-1)+1]
    newpar 2pds:[expr 3*($nlor-1)+2]=1pds:p[expr 3*($nlor-1)+2]
    newpar 2pds:[expr 3*($nlor-1)+3] $norm 0.01 1e-10 1e-10 1e10 1e10

    puts "Adding xlor $nlor to 3rea"
    addco 3rea:$nlor xlor & /*
    newpar 3rea:[expr 4*($nlor-1)+1]=1pds:p[expr 3*($nlor-1)+1]
    newpar 3rea:[expr 4*($nlor-1)+2]=1pds:p[expr 3*($nlor-1)+2]
    newpar 3rea:[expr 4*($nlor-1)+3] 0.000 0.01 -$::pi -$::pi $::pi $::pi
    newpar 3rea:[expr 4*($nlor-1)+4] $norm 0.01 1e-10 1e-10 1e10 1e10

    puts "Adding ylor $nlor to 4ima"
    addco 4ima:$nlor ylor & /*
    newpar 4ima:[expr 4*($nlor-1)+1]=3rea:p[expr 4*($nlor-1)+1]
    newpar 4ima:[expr 4*($nlor-1)+2]=3rea:p[expr 4*($nlor-1)+2]
    newpar 4ima:[expr 4*($nlor-1)+3]=3rea:p[expr 4*($nlor-1)+3]
    newpar 4ima:[expr 4*($nlor-1)+4]=3rea:p[expr 4*($nlor-1)+4]

    set lors [tcloutr modcomp 1pds]
    puts "Rearranging plag and coherence to $lors Lors"
    plags $lors
    cohes $lors
    puts "Finished adding Lor $nlor"
    chatter $chatlevel
}


# plags --
#   Creates a phase lag (plag) model for spectrum 5pla.
# Arguments:
#   args - Number of Lorentzians
proc plags {args} {
    set model_string [tcloutr model]
    if {[string match "*5pla*" $model_string]} {
        puts "Erasing existing phase lag model"
        mdefine plags :
    }
    set chatlevel [scan [tcloutr chatter] "%d"]
    chatter 0
    set nlor [lindex $args 0]
    puts "Constructing phase lag model with $nlor Lors"
    set num "("
    set den "("
    for {set i 1} {$i <= $nlor} {incr i} {
        append num "xlor(LineE$i,Width$i,plag$i)*norm$i+"
        append den "ylor(LineE$i,Width$i,plag$i)*norm$i+"
    }
    set pe "atan2([string range $den 0 end-1]),[string range $num 0 end-1])) - $::pshift"
    query yes
    mdefine plags $pe : add
    puts "Loading the 5pla plag model with $nlor Lors"
    model 5:5pla plags & /*
    for {set i 1} {$i <= 4*$nlor} {incr i} { newpar 5pla:$i =3rea:p$i }
    newpar 5pla:$i 1,-1
    puts "Finished adding 5pla model with $nlor Lors"
    chatter $chatlevel
}


# cohes --
#   Creates a coherence model for spectrum 6coh.
# Arguments:
#   args - Number of Lorentzians
proc cohes {args} {
    set model_string [tcloutr model]
    if {[string match "*6coh*" $model_string]} {
        puts "Erasing existing coherence model"
        mdefine coherence :
    }
    set chatlevel [scan [tcloutr chatter] "%d"]
    chatter 0
    set nlor [lindex $args 0]
    puts "Constructing coherence model with $nlor Lors"
    set re "("
    set im "("
    set p1 "("
    set p2 "("
    for {set i 1} {$i <= $nlor} {incr i} {
        append re "xlor(LineE$i,Width$i,plag$i)*norm$i+"
        append im "ylor(LineE$i,Width$i,plag$i)*norm$i+"
        append p1 "Lorentz(LineE$i,Width$i)*nn$i+"
        append p2 "Lorentz(LineE$i,Width$i)*mm$i+"
    }
    set pe "([string range $re 0 end-1])**2+[string range $im 0 end-1])**2)/([string range $p1 0 \
    end-1])*[string range $p2 0 end-1]))"
    query yes
    mdefine coherence $pe : add
    puts "Loading the 6coh coherence model with $nlor Lors"
    model 6:6coh coherence & /*
    for {set i 1} {$i <= 4*$nlor} {incr i} { newpar 6coh:$i = 3rea:p$i }
    for {set i 1} {$i <= $nlor} {incr i} { newpar 6coh:[expr 4*$nlor+$i] = 1pds:p[expr 3*$i] }
    for {set i 1} {$i <= $nlor} {incr i} { newpar 6coh:[expr 4*$nlor+$nlor+$i] = 2pds:p[expr 3*$i] }
    newpar 6coh:[expr 4*$nlor+2*$nlor+1] 1,-1
    puts "Finished adding 6coh coherence model with $nlor Lors"
    chatter $chatlevel
}


# sort_lorentzians --
#   Sorts Lorentzians by characteristic frequency (sqrt(freq^2 + width^2)).
proc sort_lorentzians {} {

    if {$::nlor <= 1} { return }

    set chatlevel [scan [tcloutr chatter] "%d"]
    chatter 0

    set lorentzians {}
    for {set i 1} {$i <= $::nlor} {incr i} {
        set base_idx_1pds [expr {3 * ($i - 1)}]
        set base_idx_2pds [expr {3 * ($i - 1)}]
        set base_idx_3rea [expr {4 * ($i - 1)}]
        set freq  [lindex [tcloutr param 1pds:[expr {$base_idx_1pds + 1}]] 0]
        set width [lindex [tcloutr param 1pds:[expr {$base_idx_1pds + 2}]] 0]
        set norm1 [lindex [tcloutr param 1pds:[expr {$base_idx_1pds + 3}]] 0]
        set norm2 [lindex [tcloutr param 2pds:[expr {$base_idx_2pds + 3}]] 0]
        set plag  [lindex [tcloutr param 3rea:[expr {$base_idx_3rea + 3}]] 0]
        set norm3 [lindex [tcloutr param 3rea:[expr {$base_idx_3rea + 4}]] 0]
        set char_freq [expr {sqrt($freq * $freq + $width * $width)}]
        lappend lorentzians [list $i $freq $width $norm1 $norm2 $plag $norm3 $char_freq]
    }

    set sorted_lorentzians [lsort -real -index 7 $lorentzians]
    set orig_char_freqs [lmap lor $lorentzians {lindex $lor 7}]
    set sorted_char_freqs [lmap lor $sorted_lorentzians {lindex $lor 7}]
    if {$orig_char_freqs eq $sorted_char_freqs} {
        chatter $chatlevel
        return
    }

    for {set i 0} {$i < [llength $sorted_lorentzians]} {incr i} {
        set lor [lindex $sorted_lorentzians $i]
        lassign $lor orig_idx freq width norm1 norm2 plag norm3 char_freq
        set new_base_idx_1pds [expr {3 * $i}]
        set new_base_idx_2pds [expr {3 * $i}]
        set new_base_idx_3rea [expr {4 * $i}]
        newpar 1pds:[expr {$new_base_idx_1pds + 1}] $freq
        newpar 1pds:[expr {$new_base_idx_1pds + 2}] $width
        newpar 1pds:[expr {$new_base_idx_1pds + 3}] $norm1
        newpar 2pds:[expr {$new_base_idx_2pds + 3}] $norm2
        newpar 3rea:[expr {$new_base_idx_3rea + 3}] $plag
        newpar 3rea:[expr {$new_base_idx_3rea + 4}] $norm3
    }

    chatter $chatlevel
}



# save_fit_log --
#   Saves detailed fit log with Lorentzian parameters and WW test statistics.
proc save_fit_log {basename spectra sigma} {

	set unit [open [format "%s.log" $basename] "w+"]

    puts $unit "\nDetailed Fit Log\n----------------------------------------"
    puts $unit [format "Spectra: %s" $spectra]
    puts $unit [format "Number of Lorentzians: %d" $::nlor]
    puts $unit [format "p-value limit = %e \t sigma = " $::pval $sigma]
#    puts $unit [format "Reduced chi2 threshold: %.4f" $::rchi]

    puts $unit "\nFit Statistics\n------------------------------------------"
    puts $unit [format "  Chi-squared: %f" [tcloutr stat]]
    puts $unit [format "  DOF: %d" [lindex [tcloutr dof] 0]]
    puts $unit [format "  Reduced Chi-squared: %f" \
            [expr {[tcloutr stat] / double([lindex [tcloutr dof] 0])}]]
	puts $unit  [format "  F-test:%f \t p-value:%e" $::last_ftest $::last_pvalue]

    puts $unit "\nLorentzian Parameters\n-----------------------------------"
    set headers [list "Lorentzian" "Freq" "Width" "PDS1 Norm" "PDS2 Norm" "Phase Lag" "CS Norm "]
    set widths [list 15 15 15 15 15 15 15]
    set header_line ""
    foreach header $headers width $widths {
        append header_line [format "%-*s" $width $header]
    }
    puts $unit $header_line
    puts $unit [string repeat "-" [string length $header_line]]

    for {set i 1} {$i <= $::nlor} {incr i} {
        set base_idx_ps [expr {3 * ($i - 1)}]
        set base_idx_cs [expr {4 * ($i - 1)}]
        set freq  [lindex [tcloutr param 1pds:[expr {$base_idx_ps + 1}]] 0]
        set width [lindex [tcloutr param 1pds:[expr {$base_idx_ps + 2}]] 0]
        set norm1 [lindex [tcloutr param 1pds:[expr {$base_idx_ps + 3}]] 0]
        set norm2 [lindex [tcloutr param 2pds:[expr {$base_idx_ps + 3}]] 0]
        set plag  [lindex [tcloutr param 3rea:[expr {$base_idx_cs + 3}]] 0]
        set norm3 [lindex [tcloutr param 3rea:[expr {$base_idx_cs + 4}]] 0]
        set row [list \
            [format "%d" $i] \
            [format "%.6e" $freq] \
            [format "%.6e" $width] \
            [format "%.6e" $norm1] \
            [format "%.6e" $norm2] \
            [format "%.6e" $plag] \
            [format "%.6e" $norm3] \
        ]
        set row_line ""
        foreach value $row width $widths {
            append row_line [format "%-*s" $width $value]
        }
        puts $unit $row_line
    }

    puts $unit "\nWald-Wolfowitz Runs Test Statistics (Final Fit):"
    set ww_headers [list "Spectrum" "Runs" "p-value" "Z-score" "N_pos" "N_neg"]
    set ww_widths [list 12 12 12 12 12 12]
    set ww_header_line ""
    foreach header $ww_headers width $ww_widths {
        append ww_header_line [format "%-*s" $width $header]
    }
    puts $unit $ww_header_line
    puts $unit [string repeat "-" [string length $ww_header_line]]

    if {[info exists ::last_ww_results]} {
        lassign $::last_ww_results overall_result per_spectrum_stats
        foreach stat $per_spectrum_stats {
            lassign $stat spec_name runs p z_score n_pos n_neg
            set row [list \
                $spec_name \
                [format "%d" $runs] \
                [format "%.6e" $p] \
                [format "%f" $z_score] \
                [format "%d" $n_pos] \
                [format "%d" $n_neg] \
            ]
            set row_line ""
            foreach value $row width $ww_widths {
                append row_line [format "%-*s" $width $value]
            }
            puts $unit $row_line
        }
    } else {
        puts $unit "No WW test results available."
    }

    puts $unit "\nSpectral Data:"
    set fit_data [get_fit_data $spectra]
    foreach spec $fit_data {
        lassign $spec spec_name _ residuals chi2_sum nchan
        puts $unit [format "Spectrum %s: X2=%f \t channels=%d" $spec_name $chi2_sum $nchan]
    }

	close $unit

}

