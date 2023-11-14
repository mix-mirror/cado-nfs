#!/usr/bin/env bash

old_setx=$(shopt -po xtrace || :)

set -e

exec < /dev/null

if [ "$CADO_DEBUG" ] ; then set -x ; fi

# This script is intended *for testing only*. It's used for, e.g.,
# coverage tests.  We don't try to use this script for real examples.
# Real examples require the user to craft his own bwc.pl command line
# with care (command-lines as created by this tool might be a source of
# inspiration, though).

pass_bwcpl_args=()

while [ $# -gt 0 ] ; do
    a="$1"
    shift
    if [ "$a" = "--" ] ; then
        break
    else
        if [[ $a =~ ^(prime|m|n|wdir|interval|seed|nullspace|mpi|thr|simd)= ]] ; then
            # There are quite a few parameters that we parse here and
            # pass to bwcpl anyway, so let's not put them in
            # pass_bwcpl_args.
            # (in fact, this is the case of most parameters)
            :
        # (try to remove both the exception here and the offending
        # exception later on)
        # elif [[ $a =~ ^(tolerate_failure|stop_at_step|keep_rolling_checkpoints|checkpoint_precious|skip_online_checks|interleaving) ]] ; then
        #     # basically the same story here, except that there seem to
        #     # be two ways these arguments get passed to bwc.pl : we
        #     # forcibly add them to ${common[@]} later on in this file,
        #     # but why do we do that? Is it only to provide for the case
        #     # where these arguments are exported via the shell
        #     # environment?
        #     :
        elif [[ $a =~ ^(bindir|mats|pre_wipe|random_matrix_size|random_matrix_minkernel|script_steps|nrhs|sage|wordsize) ]] ; then
            # and there are even parameters that only make sense here.
            :
        else
            pass_bwcpl_args+=("$a")
        fi
        eval "$a"
    fi
done


# various configuration variables. environment can be used to override them
: ${scriptpath=$0}
: ${m=8}
: ${n=4}
: ${prime=4148386731260605647525186547488842396461625774241327567978137}
: ${mpi:=1x1}
: ${lingen_mpi:=1x1}
: ${thr:=2x2}
# Set the "matrix" variable in order to work on a real matrix.
: ${matrix=}
# Set the "bindir" variable to use pre-built binaries (must point to the
# directories with the bwc binaries)
: ${bindir=}

# This is related to the random matrix which gets automatically created
# if no matrix was given on the cmdline.
: ${random_matrix_size=1000}
# there is no random_matrix_coeffs_per_row ; see random_matrix.c
: ${random_matrix_maxcoeff=10}
: ${random_matrix_minkernel=10}
: ${mats=$HOME/Local/mats}
: ${pre_wipe=}
: ${seed=$RANDOM}
: ${balancing_options=reorder=columns}
: ${script_steps=wipecheck,matrix,bwc.pl/:complete,bwccheck}

pass_bwcpl_args+=("seed=$seed")

: ${wordsize=64}

# XXX note that $wdir is wiped out by this script !
: ${wdir=/tmp/bwcp}

# We also have sagemath testing.
: ${sage=}

# This has to be provided as auxiliary data for the GF(p) case (or we'll
# do without any rhs whatsoever).
: ${rhs=}

# For the rest of this script, we'll prefer the variable name "rhsfile"
rhsfile="$rhs"
: ${nullspace=right}
: ${interval=50}
# : ${mm_impl=basicp}

if ! type -p seq >/dev/null ; then
    seq() {
        first="$1"
        shift
        if [ "$2" ] ; then
            incr="$1"
            shift
        else
            incr=1
        fi
        last="$1"
        while [ "$first" -le "$last" ] ; do
            echo "$first"
            let first=first+incr
        done
    }
fi

# inject the variables that were provided by guess_mpi_configs
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
    pass_bwcpl_args+=(mpi_extra_args="${mpi_extra_args[*]}")
fi

usage() {
    cat >&2 <<-EOF
    Usage: $scriptpath [param1=value param2=value2 ...]
    This runs a complete block Wiedemann algorithm, given the parameters
    specified by various environment variables (the command line can also
    be used to set those variables).
EOF
    exit 1
}

argument_checking() {
    case "$nullspace" in
        left|right) ;;
        *) echo "\$nullspace must be left or right" >&2; usage;;
    esac
    if [ "$prime" != 2 ] ; then
        if [ "$matrix" ] && [ "$rhsfile" ] && [ "$nrhs" ] ; then
            # inhomogeneous mod p
            :
        elif [ "$matrix" ] && ! [ "$rhsfile" ] && ! [ "$nrhs" ] ; then
            # homogeneous mod p (bogus as of 755e9c5).
            :
        elif ! [ "$matrix" ] && ! [ "$rhsfile" ] ; then
            if ! [ "$random_matrix_size" ] || ! [ "$random_matrix_maxcoeff" ] ; then
                echo "Please give random matrix dimensions" >&2
                exit 1
            fi
            if ! [ "$nrhs" ] && ! [ "$random_matrix_minkernel" ] ; then
                echo "Please give random matrix dimensions" >&2
                exit 1
            fi
        else
            echo "Unsupported combination of arguments for specifying system" >&2
            echo "Detected arguments:" >&2
            for a in matrix rhsfile nrhs random_matrix_size random_matrix_maxcoeff ; do
                if [ "${!a}" ] ; then echo "$a=${!a}" >&2 ; fi
            done
            echo "Supported combinations:">&2
            echo "matrix rhsfile nrhs" >&2
            echo "matrix" >&2
            echo "random_matrix_size random_matrix_maxcoeff random_matrix_minkernel" >&2
            echo "random_matrix_size random_matrix_maxcoeff nrhs" >&2
            usage
        fi
    fi
}

derived_variables() {
    if ! [ -d $mats ] ; then mats=$wdir; fi
    if [ "$prime" = 2 ] ; then
        splitwidth=64
    else
        splitwidth=1
    fi
    : ${simd=$splitwidth}
}

prepare_wdir() {
    if [[ $script_steps =~ wipecheck ]] ; then
        if [ -d $wdir ] ; then
            if [ "$pre_wipe" ] ; then
                rm -rf $wdir 2>/dev/null
            else
                echo "Won't wipe $wdir unless \$pre_wipe is set" >&2
                exit 1
            fi
        fi
        mkdir $wdir
    elif [[ $script_steps =~ keepdir ]] ; then
        if ! [ -d $wdir ] ; then
            echo "cannot work with keepdir: $wdir does not exist" >&2
            exit 1
        fi
    else
        echo "Want either wipecheck or keepdir in script_steps" >&2
        exit 1
    fi
}


create_test_matrix_if_needed() {
    if [ "$matrix" ] ; then
        if ! [ -e "$matrix" ] ; then
            echo "matrix=$matrix inaccessible" >&2
            usage
        fi
        # Get absolute path.
        matrix=$(readlink /proc/self/fd/99 99< $matrix)
        if [ -e "$rhsfile" ] ; then
            rhsfile=$(readlink /proc/self/fd/99 99< $rhsfile)
        fi
        return
    fi

    # It's better to look for a kernel which is not trivial. Thus
    # specifying --kright for random generation is a good move prior to
    # running this script for nullspace=right
    kside="--k$nullspace"

    # We only care the binary matrix, really. Nevertheless, the random
    # matrix is created as text, and later transformed to binary format.
    # We also create the auxiliary .rw and .cw files.

    # We're not setting the density, as there is an automatic setting in
    # random_matrix for that (not mandatory though, since
    # nrows,ncols,density, may be specified in full).

    rmargs=()
    # defaults, some of the subcases below tweak that.
    nrows=`echo $random_matrix_size | cut -d, -f1`
    ncols=`echo $random_matrix_size | cut -d, -f2`   # = nrows if no comma
    outer_nrows=`echo $random_matrix_size | cut -d, -f1`
    outer_ncols=`echo $random_matrix_size | cut -d, -f2`   # = nrows if no comma
    escaped_size=$(echo $random_matrix_size | tr , _)
    if [ "$prime" = 2 ] ; then
        basename=$mats/t${escaped_size}
        matrix="$basename.matrix.bin"
        rmargs+=(--k$nullspace ${random_matrix_minkernel})
        # ncols=
    elif ! [ "$nrhs" ] ; then
        basename=$mats/t${escaped_size}p
        matrix="$basename.matrix.bin"
        rmargs+=(--k$nullspace ${random_matrix_minkernel})
        rmargs+=(-c ${random_matrix_maxcoeff})
        rmargs+=(-Z)
    else
        # This is an experimental mode. In the DLP context, we have a
        # matrix with N rows, N-r ideal columns, and r Schirokauer maps.   
        # let's say random_matrix_size is N and nrhs is r. We'll generate
        # a matrix with N rows and N-r columns, and later pad the column
        # width data with r zeroes.
        ncols=$((nrows-nrhs))
        basename=$mats/t${escaped_size}p+${nrhs}
        # We want something proportional to log(N)^2 as a density.
        # There's no undebatable support for this heuristic, it's just an
        # arbitrary choice. We'll do that without relying on external
        # tools for computing the log: just plain stupid shell. That's
        # gonna be good enough. The ratio below is taken as 2, which fits
        # well the relevant density ranges for the matrices we wish to
        # consider.
        density=$((${#random_matrix_size}*${#random_matrix_size}*2))
        if [ "$density" -lt 12 ] ; then density=12; fi
        rmargs+=(-d $density -Z)
        matrix="$basename.matrix.bin"
        rhsfile="$basename.rhs.txt"
        rmargs+=(-c ${random_matrix_maxcoeff})
        rmargs+=(rhs="$nrhs,$prime,$rhsfile")
    fi
    rmargs=($nrows $ncols -s $seed "${rmargs[@]}" --freq --binary --output "$matrix")
    rwfile=${matrix%%bin}rw.bin
    cwfile=${matrix%%bin}cw.bin
    ncols=$outer_ncols
    nrows=$outer_nrows
    if [[ $script_steps =~ matrix ]] ; then
        ${bindir}/random_matrix "${rmargs[@]}"
        data_ncols=$((`wc -c < $cwfile` / 4))
        data_nrows=$((`wc -c < $rwfile` / 4))
        if [ "$data_ncols" -lt "$ncols" ] ; then
            if [ "$prime" = 2 ] || ! [ "$nrhs" ] ; then
                echo "padding $cwfile with $((ncols-data_ncols)) zero columns"
                dd if=/dev/zero bs=4 count=$((ncols-data_ncols)) >> $cwfile
            fi
        fi
        if [ "$data_nrows" -lt "$nrows" ] ; then
            echo "padding $cwfile with $((nrows-data_nrows)) zero rows"
            dd if=/dev/zero bs=4 count=$((nrows-data_nrows)) >> $rwfile
        fi
    fi
}

# This is only useful in the situation where an input matrix has been
# provided by the user (which may be ither in text or binary format). In
# this case, we must make sure that we have the .cw and .rw files too.
create_auxiliary_weight_files() {
    if ! [[ $script_steps =~ matrix ]] ; then
        return
    fi
    if [ "$prime" != 2 ] ; then withcoeffs=--withcoeffs ; fi
    case "$matrix" in
        *.txt)
            if [ "$rhsfile" ] ; then
                # It's really a hassle to keep the conversion code (which
                # existed until dad7019)
                echo "Please supply $rhs as a *binary* file, please\n" >&2
                exit 1
            fi
            matrix_txt="$matrix"
            matrix=${matrix%%txt}bin
            rwfile=${matrix%%bin}rw.bin
            cwfile=${matrix%%bin}cw.bin
            if [ "$matrix" -nt "$matrix_txt" ] && [ "$rwfile" -nt "$matrix_txt" ] && [ "$cwfile" -nt "$matrix_txt" ] ; then
                echo "Taking existing $mfile, $rwfile, $cwfile as accompanying $matrix_txt"
            else
                echo "Creating files $matrix, $rwfile, $cwfile from $matrix_txt"
                $bindir/mf_scan  --ascii-in $withcoeffs --mfile $matrix_txt  --freq --binary-out --ofile $matrix
            fi
            ;;
        *.bin)
            rwfile=${matrix%%bin}rw.bin
            cwfile=${matrix%%bin}cw.bin
            if [ "$rwfile" -nt "$matrix" ] && [ "$cwfile" -nt "$matrix" ] ; then
                echo "Taking existing $rwfile, $cwfile as accompanying $matrix"
            else
                $bindir/mf_scan  --binary-in $withcoeffs --mfile $matrix --freq
            fi
            ;;
    esac
    ncols=$((`wc -c < $cwfile` / 4))
    nrows=$((`wc -c < $rwfile` / 4))
}

prepare_common_arguments() {
    # This sets the common arguments for all bwc binaries
    common=(
        matrix=$matrix
        balancing_options=$balancing_options
        mpi=$mpi
        thr=$thr
        m=$m
        n=$n
        wdir=$wdir
        prime=$prime
        nullspace=$nullspace
        interval=$interval
        simd=$simd
    )
    if [ "$mm_impl" ] ; then common+=(mm_impl=$mm_impl) ; fi
    if [ "$prime" != 2 ] ; then common+=(lingen_mpi=$lingen_mpi) ; fi
}

argument_checking
derived_variables
prepare_wdir

if ! [ "$matrix" ] ; then
    # This also sets rwfile cwfile nrows ncols
    create_test_matrix_if_needed
else
    # This also sets rwfile cwfile nrows ncols
    create_auxiliary_weight_files
fi
for f in matrix cwfile rwfile ; do
    if ! [ -f "${!f}" ] ; then echo "Missing file $f=${!f}" >&2 ; exit 1 ; fi
done

prepare_common_arguments

if [ "$rhsfile" ] ; then
    common+=(rhs=$rhsfile)
fi

# for v in tolerate_failure stop_at_step keep_rolling_checkpoints checkpoint_precious skip_online_checks interleaving ; do
#     if [ "${!v}" ] ; then common+=("$v=${!v}") ; fi
# done

if [ "$sage" ] ; then
    common+=(save_submatrices=1)
fi

if [ "$sage" ] ; then
    echo "### Enabling SageMath checking ###"
fi

if ! [ "$sage" ] ; then
    rc=0
    if [[ $script_steps =~ bwc\.pl/([a-z:/]*) ]] ; then
        IFS=/ read -a bwcpl_steps <<< "${BASH_REMATCH[1]}"
        for s in "${bwcpl_steps[@]}" ; do
            if $bindir/bwc.pl "$s" "${common[@]}" "${pass_bwcpl_args[@]}" ; then
                :
            else
                rc=$?
                break
            fi
        done
    fi
    if [[ $script_steps =~ bwccheck ]] ; then
        $bindir/bwc.pl :mpirun_single -- $bindir/bwccheck prime=$prime m=$m n=$n -- $wdir/[ACVFS]* > $wdir/bwccheck.log
        grep NOK $wdir/bwccheck.log || :
    fi
    exit $rc
else
    set +e
    rc=0
    if [[ $script_steps =~ bwc\.pl/([a-z:/]*) ]] ; then
        IFS=/ read -a bwcpl_steps <<< "${BASH_REMATCH[1]}"
        for s in "${bwcpl_steps}" ; do
            if $bindir/bwc.pl "$s" "${common[@]}" "${pass_bwcpl_args[@]}" ; then
                :
            else
                rc=$?
                break
            fi
        done
    fi
    if [[ $script_steps =~ bwccheck ]] ; then
        $bindir/bwc.pl :mpirun_single -- $bindir/bwccheck prime=$prime m=$m n=$n -- $wdir/[ACVFS]* > $wdir/bwccheck.log
        grep NOK $wdir/bwccheck.log
    fi
    eval $old_setx
    if [[ $script_steps =~ bwc\.pl([a-z:/]*) ]] ; then
        if [ "$rc" = 0 ] ; then
            echo " ========== SUCCESS ! bwc.pl returned true ========== "
            echo " ========== SUCCESS ! bwc.pl returned true ========== "
            echo " ========== SUCCESS ! bwc.pl returned true ========== "
        else
            echo " ########## FAILURE ! bwc.pl returned false ########## "
            echo " ########## FAILURE ! bwc.pl returned false ########## "
            echo " ########## FAILURE ! bwc.pl returned false ########## "
        fi
    fi
fi

eval $old_setx

sage_check_parameters() { # {{{
    # Required parameters below this point:

    # wdir
    # matrix
    # m
    # n
    # prime
    # interval
    # splitwidth
    # cmd (pay attention to dirname $0 above !)
    # rwfile
    # cwfile
    # nullspace

    for v in wdir matrix m n prime interval splitwidth cmd rwfile cwfile nullspace ; do
        if ! [ "${!v}" ] ; then echo "Missing parameter \$$v" >&2 ; fi
    done

    # nrows and ncols are also needed, but can be deduced from rwfile and
    # cwfile easily.
    : ${ncols=$((`wc -c < $cwfile` / 4))}
    : ${nrows=$((`wc -c < $rwfile` / 4))}

    mdir=$wdir

    if ! [[ "$mpi,$thr" =~ ^([0-9]*)x([0-9]*),([0-9]*)x([0-9]*)$ ]] ; then
        echo "bad format for mpi=$mpi and thr=$thr" >&2
        exit 1
    fi

    Nh=$((${BASH_REMATCH[1]}*${BASH_REMATCH[3]}))
    Nv=$((${BASH_REMATCH[2]}*${BASH_REMATCH[4]}))

    bfile="`basename $matrix .bin`.${Nh}x${Nv}/`basename $matrix .bin`.${Nh}x${Nv}.bin"

    : ${nrhs:=0}
}
# }}}

if [ "$sage" ] ; then
    cmd=/bin/true
    sage_check_parameters
    cd `dirname "$0"`
    export PYTHONUNBUFFERED=true
    check_script_diagnostic_fd=1
    if [ "$FORCE_BWC_EXTERNAL_CHECKS_OUTPUT_ON_FD3" ] && (exec 1>&3) 2>&- ; then
        check_script_diagnostic_fd=3
    fi
    sage_args=( m=$m n=$n p=$prime
                wdir=$wdir matrix=$matrix
                nh=$Nh nv=$Nv
    )
    if [ "$wordsize" != 64 ] ; then
        sage_args+=(wordsize=$wordsize)
    fi
    if [ "$CADO_DEBUG" ] ; then set -x ; fi
    set -eo pipefail
    "$sage" bwc.sage "${sage_args[@]}" >&${check_script_diagnostic_fd}
    eval $old_setx
fi
