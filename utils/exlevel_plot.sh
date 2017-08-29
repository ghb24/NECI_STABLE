#!/bin/bash
#-------------------------------------------------------#
# Plot data from EXLEVELStats as either an animated bar #
# graph or a single lineplot.                           #
#                                                       #
# PLR 03.2017                                           #
#-------------------------------------------------------#
set +u
shopt -s extglob


########################### UTILITIES ###########################


repeat() {
  #--------------------#
  # Print $2 $1 times. #
  #--------------------#
  local i string
  string=""
  for ((i=1; i<=$1; i++)) ; do
    string="$string$2"
  done
  echo "$string"
}


zero_pad() {
  #----------------------------------------------------#
  # Pad integer $2 with zeroes until it has length $1. #
  #----------------------------------------------------#
  echo "$(repeat $(($1-${#2})) 0)$2"
}


errstop() {
  #-----------------------------#
  # Stop with error message $1. #
  #-----------------------------#
  echo "$EL"
  echo "ERROR: $1"
  echo
  exit 1
}


check_number_N() {
  #-------------------------------------------------------#
  # Return with status 0 if $1 is a non-negative integer. #
  #-------------------------------------------------------#
  [[ "$1" == +([[:digit:]]) ]]
}


check_number_R() {
  #----------------------------------------------#
  # Return with status 0 if $1 is a real number. #
  #----------------------------------------------#
  [[ "$1" ==\
     ?([+-])+([[:digit:]])?(.)*([[:digit:]])?([eE]?([+-])+([[:digit:]])) ]]\
     || [[ "$1" ==\
     ?([+-])*([[:digit:]]).+([[:digit:]])?([eE]?([+-])+([[:digit:]])) ]]
}


float2human() {
  #----------------------------------------#
  # Format number $1 in a readable manner. #
  #----------------------------------------#
  local n i p
  n=$(printf "%f" "$1")
  n=${n%.*}
  ((i=n))
  p=0
  if ((n<1000)) ; then
    o=$n
  else
    while ((i>=1000)) ; do
      ((i/=1000))
      ((p+=3))
    done
    if ((i>99)) ; then
      o=${n:0:3}
    elif ((i>9)) ; then
      o=${n:0:2}.${n:2:1}
    else
      o=${n:0:1}.${n:1:2}
    fi
  fi
  case "$p" in
  0) echo "$o" ;;
  3) echo "${o}k" ;;
  6) echo "${o}M" ;;
  9) echo "${o}G" ;;
  12) echo "${o}T" ;;
  15) echo "${o}P" ;;
  18) echo "${o}E" ;;
  21) echo "${o}Z" ;;
  24) echo "${o}Y" ;;
  *) echo "${o}+$p" ;;
  esac
}


########################### FUNCTIONS ###########################


general_setup() {
  #-----------------#
  # General set-up. #
  #-----------------#

  # Get tput codes for delete-to-end-of-line and carriage return.
  [[ "$TERM" == xterm-* ]] && export TERM=xterm
  EL=$(tput el 2> /dev/null)
  [ -z "$EL" ] && el=$(tput ce 2> /dev/null)
  CR=$(tput cr 2> /dev/null)

}


read_command_line() {
  #-----------------------------------------#
  # Set options from the command line "$@". #
  #-----------------------------------------#
  local option option1 var val minval maxval
  local delete_framedir gnuplot_version image_width video_format

  # Print header.
  echo "EXLEVEL PLOTTING TOOL"
  echo "====================="

  # Set defaults.
  TYPE=bar-video
  RELATIVE=1
  REF_LAST=0
  delete_framedir=0
  EXSTATS=EXLEVELStats
  STATS=FCIMCStats
  FRAMEDIR=./frames
  REF_DATA=""
  GNUPLOT=gnuplot
  FFMPEG=ffmpeg
  ! type -P ffmpeg >& /dev/null && type -P avconv >& /dev/null && FFMPEG=avconv
  video_format=avi
  LNORM=1
  NFRAME=0
  image_width=640
  FRAME_RATE=10
  DRAW_EVERY=1

  # Parse command line.
  while (($#>0)) ; do

    case "$1" in

    --*) # GNU-style "long" options
      option="${1#--}"
      case "$option" in
      absolute)
        RELATIVE=0 ;;
      ref-last)
        REF_LAST=1 ;;
      force)
        delete_framedir=1 ;;
      gif)
        video_format=gif ;;
      plot)
        TYPE=line-plot ;;
      video)
        TYPE=bar-video ;;
      frames)
        TYPE=bar-frames ;;
      help)
        dump_help
        exit ;;

      *=*) # value-setting long options
        var="${option%%=*}" ; val="${option#*=}"

        case "$var" in
        exstats|stats|framedir|ref-data|gnuplot|ffmpeg) # strings
          # Regularize variable name.
          var="${var^^}"
          var="${var//-/_}"
          # Set value.
          eval "$var=\"\$val\"" ;;

        lnorm|nframe|width|frame-rate|draw-every) # integers
          check_number_N $val || errstop\
             "Argument to --$var must be an integer."
          # Renames.
          case "$var" in
          width) var=image_width ;;
          *)
            # Regularize variable name.
            var="${var^^}"
            var="${var//-/_}" ;;
          esac
          # Set bounds.
          minval=""
          maxval=""
          case "$var" in
          LNORM)       minval=0 ; maxval=2 ;;
          NFRAME)      minval=1 ;;
          image_width) minval=1 ;;
          FRAME_RATE)  minval=1 ;;
          DRAW_EVERY)  minval=1 ;;
          esac
          # Check bounds.
          [ ! -z "$minval" ] && ((val<minval)) && errstop\
             "Argument to --$var must be greater than or equal to $minval."
          [ ! -z "$maxval" ] && ((val>maxval)) && errstop\
             "Argument to --$var must be less than or equal to $maxval."
          # Set value.
          eval "$var=\"\$val\"" ;;

        *) errstop "Unrecognized option --$option." ;;
        esac ;;

      *) errstop "Unrecognized option --$option."
      esac ;;

    -*) # Unix-style "short" options
      option=${1#-}
      [ -z "$option" ] && errstop "Bad option '-'."

      while ((${#option}>0)) ; do
        option1=${option:0:1}
        option=${option:1}

        case "$option1" in
        a)
          RELATIVE=0 ;;
        L)
          REF_LAST=1 ;;
        f)
          delete_framedir=1 ;;
        p)
          TYPE=line-plot ;;
        V)
          TYPE=bar-video ;;
        v)
          TYPE=bar-frames ;;
        h)
          dump_help
          exit ;;

        X|S|D|R|G|F) # strings
          if [ ! -z "$option" ] ; then
            val="$option"
            option=""
          else
            (($#==1)) && errstop "-$option1 must be followed by a string."
            shift
            val="$1"
          fi
          # Set variable names.
          case "$option1" in
          X) var=EXSTATS ;;
          S) var=STATS ;;
          D) var=FRAMEDIR ;;
          R) var=REF_DATA ;;
          G) var=GNUPLOT ;;
          F) var=FFMPEG ;;
          esac
          # Set variable.
          eval "$var=\"\$val\"" ;;

        l|n|w|r|d) # integers
          if [ ! -z "$option" ] ; then
            val="$option"
            option=""
          else
            (($#==1)) && errstop "-$option1 must be followed by an integer."
            shift
            val="$1"
          fi
          # Check value type.
          check_number_N "$val" || errstop\
             "Argument to -$option1 must be an integer."
          # Set variable names and value bounds.
          minval=""
          maxval=""
          case "$option1" in
          l) var=LNORM       ; minval=0 ; maxval=2 ;;
          n) var=NFRAME      ; minval=1 ;;
          w) var=image_width ; minval=1 ;;
          r) var=FRAME_RATE  ; minval=1 ;;
          d) var=DRAW_EVERY  ; minval=1 ;;
          esac
          # Check bounds.
          [ ! -z "$minval" ] && ((val<minval)) && errstop \
             "Argument to -$option1 must be greater than or equal to $minval."
          [ ! -z "$maxval" ] && ((val>maxval)) && errstop \
             "Argument to -$option1 must be less than or equal to $maxval."
          # Set variable.
          eval "$var=\"\$val\"" ;;

        *) errstop "Unrecognized option -$option1." ;;
        esac

      done ;;

    *)
      errstop "Unrecognized argument $option" ;;
    esac
    shift

  done

  # Check for required files.
  [ -e "$STATS" ] || errstop "File '$STATS' not found."
  [ -e "$EXSTATS" ] || errstop "File '$EXSTATS' not found."
  [ ! -z "$REF_DATA" ] && [ ! -e "$REF_DATA" ]\
     && errstop "File '$REF_DATA' not found."

  # Check for required programs.
  type -P "$GNUPLOT" >& /dev/null || errstop "'$GNUPLOT' binary not found."
  gnuplot_version=$($GNUPLOT --version 2> /dev/null)
  gnuplot_version="$(set -- $gnuplot_version ; echo $2)"
  case "$gnuplot_version" in
  1.*|2.*|3.*|4.0*|4.1*)
    errstop "Gnuplot version $gnuplot_version not supported." ;;
  esac
  if [ "$TYPE" = bar-video ]  ; then
    type -P "$FFMPEG" >& /dev/null || errstop "'$FFMPEG' binary not found."
  fi

  # Set filename stem.
  STEM=EXLEVEL_W_L$LNORM
  case "$TYPE" in
  bar-video|bar-frames)
    # This is a multi-image type.
    VIDEO=$STEM.$video_format
    # Set resolution.
    ((X_RESOLUTION=image_width))
    ((Y_RESOLUTION=(3*image_width)/4)) # enforce 4:3 aspect ratio
    # Create FRAMEDIR.
    FRAMEDIR="$FRAMEDIR/L$LNORM"
    if [ -d "$FRAMEDIR" ] ; then
      if ((delete_framedir)) ; then
        rm -rf "$FRAMEDIR" >& /dev/null
      else
        errstop "Directory '$FRAMEDIR' already present.  Use '-f' to remove."
      fi
    fi
    mkdir -p "$FRAMEDIR" >& /dev/null ;;
  line-plot)
    # This is a single-image type.
    PLOT=$STEM.gpi ;;
  esac

}


dump_help() {
  #-------------------------#
  # Print help information. #
  #-------------------------#

  cat <<__EOF
Plot data from EXLEVELStats as either an animated bar graph or a single plot.
Requires gnuplot 5+, and either ffmpeg or avconv for animated bar graphs.

Usage
-----
exlevel_plot.sh [<optional-arguments>]

Command-line arguments
----------------------
--bar-video | -V   [default]
--bar-frames | -v
--line-plot | -p
  These options set the output type to be an animated bar graph, the frames
  in the animated bar graph, or a single line plot, respectively.

--absolute | -a
  By default this script plots L-norms of the weights of each excitation
  level relative to that of the zeroeth excitation level.  This option
  causes absolute L-norms of the weights to be plotted instead.

--lnorm=<l-norm> | -l <l-norm>
  Requests plotting L<l-norm> norm of the weights.  L0, L1 and L2 norms
  are available in the EXLEVELStats file.
  Default: 1

--exstats=<exstats> | -X <exstats>
  Sets name of the EXLEVELStats file.
  Default: EXLEVELStats

--stats=<stats> | -S <stats>
  Sets name of the FCIMCStats file.
  Default: FCIMCStats

--gnuplot=<gnuplot> | -G <gnuplot>
  Sets the name of the gnuplot binary.  Note that gnuplot is used during data
  processing to determine plot ranges, so it is mandatory even in line-plot
  mode.
  Default: gnuplot

--help | -h
  Shows this help.

Command-line arguments for bar-video and bar-frames modes
---------------------------------------------------------
--ref-data=<ref-data> | -R <ref-data>
  Plots reference data file <ref-data> in every frame of the animation.
  This should be a .dat file under a frames/ directory produced by a
  previous execution of this script.
  Default: (none)

--ref-last | -L
  Uses the data for the last frame of the animation as reference.  Overrides
  the use of --ref-data.

--width=<width> | -w <width>
  Sets the width in pixels of the video frames.  NB, the aspect ratio is
  fixed at 4:3.
  Default: 640

--framedir=<frame-dir> | -D <frame-dir>
  Sets directory for storing frame data and PNGs.
  Default: ./frames

--force | -f
  This script aborts with an error if the frames directory is already
  present.  This option forces the deletion of exisiting frames.

Command-line arguments for bar-video mode
-----------------------------------------
--draw-every=<draw-every> | -d <draw-every>
  Skip <draw-every>-1 data lines per frame, for shorter videos.
  Default: 1

--nframe=<nframe> | -n <nframe>
  Trim the animation to <nframe> frames.
  Default: the number of data lines in EXLEVELStats divided by <draw-every>

--frame-rate=<frame-rate> | -r <frame-rate>
  Sets the video frame rate in frames per second.
  Default: 10

--ffmpeg=<ffmpeg> | -F <ffmpeg>
  Sets the name of the ffmpeg binary.
  Default: whichever of "ffmpeg" or "avconv" is in the user's PATH.

--gif
  Produce a .gif file instead of a .avi file.  GIF players tend to offer no
  playback controls, so AVI files are generated by default.

__EOF

}


parse_files() {
  #-----------------------------------------------------------#
  # Get NEL, NFRAME, SHIFT_ACTIVE, STEP_WALKERS, STEP_ENERGY, #
  # STEP_TAU, PLOT_MINX, PLOT_MAXX, PLOT_MINY, and PLOT_MAXY. #
  #-----------------------------------------------------------#
  local line iframe istep s n t e dum miny icol0 icol iex plot_u nstep active
  local y1 y2 label val

  # Report.
  echo "Set up:$EL"
  echo "* L${LNORM} norm requested.$EL"

  # Get NEL from STATS file.
  NEL=0
  read line < $STATS
  set -- $line
  while (($#>0)) ; do
    case "$1" in
    NEl=)
      shift
      NEL=$1
      break ;;
    NEl=*)
      NEL=${1#NEl=}
      break ;;
    esac
    shift
  done
  echo "* Excitation levels: 0 -- $NEL$EL"

  # Get NFRAME from EXSTATS file.
  ((NFRAME<1)) && ((NFRAME=$(grep -cvE "^ *#" "$EXSTATS")/DRAW_EVERY))

  # Get SHIFT_ACTIVE and nstep from STATS file.
  echo -n "Parsing $STATS...$EL$CR"
  ((nstep=0))
  SHIFT_NREGION=0
  ((active=0))
  while read istep s dum dum n dum dum dum dum dum dum dum dum dum dum dum\
     dum dum t dum dum dum e dum; do
    [ "${istep:0:1}" = '#' ] && continue
    ((istep>nstep)) && ((nstep=istep))
    if [ "$s" = "0.000000" ] ; then
      if ((active)) ; then
        ((SHIFT_REGION_END[SHIFT_NREGION]=istep))
        ((active=0))
      fi
    else
      if ! ((active)) ; then
        ((SHIFT_NREGION++))
        ((SHIFT_REGION_START[SHIFT_NREGION]=istep))
        ((active=1))
      fi
      ((SHIFT_ACTIVE[istep]=1))
    fi
    STEP_WALKERS[$istep]="$(float2human $n)"
    STEP_ENERGY[$istep]="$e"
    STEP_TAU[$istep]="$t"
  done < "$STATS"
  ((active)) && ((SHIFT_REGION_END[SHIFT_NREGION]=nstep))
  echo "* $nstep MC steps in calculation, $NFRAME data points to plot.$EL"
  echo "* Number of times shift is active: $SHIFT_NREGION$EL"

  # Find overall maximum and non-zero minimum in data.
  echo -n "Parsing $EXSTATS...$EL$CR"
  miny=""
  maxy=""
  ((icol0=2+LNORM*(NEL+1)))
  for ((iex=0; iex<=NEL; iex++)) ; do
    ((icol=icol0+iex))
    plot_u="1:(\$$icol>0 ? \$$icol : 1/0)"
    ((RELATIVE)) && plot_u="1:(\$$icol>0 ? \$$icol/\$$icol0 : 1/0)"
    y1=""
    y2=""
    {
      while read label val ; do
        case "$label" in
        Y_MIN)
          y1=$val ;;
        Y_MAX)
          y2=$val ;;
        esac
      done
    } < <( $GNUPLOT 2>&1 <<____EOF
set term unknown
if ( GPVAL_VERSION == "4.2" || GPVAL_VERSION == "4.3" ) \
     set logscale y ;\
     plot '$EXSTATS' u $plot_u w l ;\
     print "Y_MIN ", GPVAL_Y_MIN ;\
     print "Y_MAX ", GPVAL_Y_MAX ;\
   else \
     plot '$EXSTATS' u $plot_u w l ;\
     print "Y_MIN ", GPVAL_DATA_Y_MIN ;\
     print "Y_MAX ", GPVAL_DATA_Y_MAX
____EOF
    )
    (exit 0) # work around fd leak in buggy versions of bash
    # Skip if y1 is not a valid number.
    [ ! -z "$y1" ] && check_number_R "$y1" && [ "$y1" != "0.0" ] || continue
    EX_PRESENT[$iex]=1
    if [ -z "$miny" ] ; then
      miny=$y1
      maxy=$y2
    else
      y1=$(printf "%20.20f" "$y1")
      y2=$(printf "%20.20f" "$y2")
      set -- $(bc -l <<<"if($miny<$y1) $miny else $y1 ;\
                         if($maxy>$y2) $maxy else $y2")
      miny=$1
      maxy=$2
    fi
  done
  echo "* Data range: $miny -- $maxy$EL"

  # Set x and y plot ranges.
  case "$TYPE" in
  bar-*)
    PLOT_MINX=-0.5
    PLOT_MAXX=$NEL.5 ;;
  line-plot)
    PLOT_MINX=0
    PLOT_MAXX=$nstep ;;
  esac
  PLOT_MINY=$(bc -l <<<"0.5*$miny")
  PLOT_MAXY=$(bc -l <<<"2*$maxy")

}


extract_frame_data() {
  #-----------------------------------------#
  # Extract frame data to individual files, #
  # and set NAME_PADDING and FRAME_STEP[:]. #
  #-----------------------------------------#
  local istep iframe jframe ilevel fbase iex y0 fdat

  # Report progress.
  echo -n "Extracting frame data...$EL$CR"

  # Loop over lines in EXSTATS file.
  NAME_PADDING=${#NFRAME}
  iframe=0
  jframe=0
  while read line ; do
    set -- $line
    istep=$1
    [ "$istep" = '#' ] && continue
    ((istep==0)) && continue
    ((jframe++))
    ((jframe%DRAW_EVERY==0)) || continue
    ((iframe++))
    echo -n "Extracting frame data... $iframe/$NFRAME$EL$CR"
    ((FRAME_STEP[iframe]=istep))
    shift $((1+LNORM*(NEL+1)))
    fdat="$FRAMEDIR/${STEM}_$(zero_pad $NAME_PADDING $iframe).dat"
    y0=$1
    for ((iex=0; iex<=NEL; iex++)) ; do
      echo "$iex $1 $y0"
      shift
    done > "$fdat"
    ((iframe>=NFRAME)) && break
  done < "$EXSTATS"

  # Adjust NFRAME downwards if needed.
  ((iframe<NFRAME)) && ((NFRAME=iframe))

  # Set reference data to last frame if so requested.
  ((REF_LAST)) && REF_DATA="$FRAMEDIR/${STEM}_$(zero_pad $NAME_PADDING\
     $NFRAME).dat"

  # Report progress.
  echo "Extracted data for $NFRAME frames.$EL"

}


gen_bar_frames() {
  #------------------#
  # Generate frames. #
  #------------------#
  local labelpos plot_u plot_ref iframe f fdat fpng
  local xlabel ylabel label_right label_left1 label_left2 istep

  # Report progress.
  echo -n "Generating frames...$EL$CR"

  # Prepare common strings to be used in gnuplot.
  labelpos="-0.2,$(bc -l <<<"0.95*$PLOT_MAXY")"
  xlabel="Excitation level"
  ylabel="L$LNORM norm of weights"
  ((RELATIVE)) && ylabel="Relative $ylabel"
  plot_u="1:2"
  plot_ref=""
  if ((RELATIVE)) ; then
    plot_u="1:(\$2/\$3)"
    [ ! -z "$REF_DATA" ] &&\
       plot_ref="'$REF_DATA' u (\$1+0.25):(\$2/\$3) w boxes ti '',"
  fi

  # Loop over frames.
  for ((iframe=1; iframe<=NFRAME; iframe++)) ; do
    # Report progress.
    echo -n "Generating frames... $iframe/$NFRAME$EL$CR"
    # Build frame.
    f="$FRAMEDIR/${STEM}_$(zero_pad $NAME_PADDING $iframe)"
    fdat="$f.dat"
    fpng="$f.png"
    ((istep=FRAME_STEP[iframe]))
    label_left1="Step #$istep | t = ${STEP_TAU[$istep\
       ]} | ${STEP_WALKERS[$istep]} walkers"
    label_left2="E = ${STEP_ENERGY[$istep]} a.u."
    label_right=""
    ((SHIFT_ACTIVE[istep])) && label_right="Shift active"
    $GNUPLOT <<____EOF
set terminal png size $X_RESOLUTION,$Y_RESOLUTION
set output '$fpng'
set boxwidth 0.4
set style fill solid
set label "$label_left1" at graph 0.05, graph 0.95 left
set label "$label_left2" at graph 0.05, graph 0.90 left
set label "$label_right" at graph 0.95, graph 0.95 right
set logscale y
set xlabel "$xlabel"
set ylabel "$ylabel"
plot [$PLOT_MINX:$PLOT_MAXX][$PLOT_MINY:$PLOT_MAXY]\
    $plot_ref '$fdat' u $plot_u w boxes ti ''
____EOF
  done

  # Report progress.
  echo "Generated $NFRAME frames.$EL"

}


gen_bar_video() {
  #-----------------#
  # Generate video. #
  #-----------------#

  # Do not generate a video file with a single frame.
  if ((NFRAME==1)) ; then
    echo "Single-frame animation, skipping video generation.$EL"
    return
  fi

  # Generate video.
  echo -n "Generating video file $VIDEO...$EL$CR"
  rm -f "$VIDEO" >& /dev/null
  if $FFMPEG -r $FRAME_RATE -i "$FRAMEDIR/${STEM}"_%0${NAME_PADDING}d.png\
     -q:v 1 -r $FRAME_RATE -s ${X_RESOLUTION}x$Y_RESOLUTION "$VIDEO" \
     >& /dev/null ; then
    echo "Video file $VIDEO generated.$EL"
  else
    echo "Problem generating $VIDEO.$EL"
  fi

}


gen_plot() {
  #--------------------------------------------------#
  # Generate gnuplot input to plot the EXSTATS file. #
  #--------------------------------------------------#
  local iregion icol0 iex icol plot_u plot_string xlabel ylabel

  # Initialize.
  echo -n "Generating $PLOT...$EL$CR"
  rm -f "$PLOT" >& /dev/null
  touch "$PLOT"
  echo "set logscale y" >> "$PLOT"

  # Draw rectangles to identify regions of variable shift.
  for ((iregion=1; iregion<=SHIFT_NREGION; iregion++)) ; do
    echo "set object $iregion rect\
       from ${SHIFT_REGION_START[$iregion]},$PLOT_MINY\
       to ${SHIFT_REGION_END[$iregion]},$PLOT_MAXY\
       fc rgb 'grey90' fs solid 1 noborder" >> "$PLOT"
  done

  # Construct plot command.
  ((icol0=2+LNORM*(NEL+1)))
  plot_string="plot [$PLOT_MINX:$PLOT_MAXX][$PLOT_MINY:$PLOT_MAXY]"
  for ((iex=0; iex<=NEL; iex++)) ; do
    ((EX_PRESENT[iex])) || continue
    ((icol=icol0+iex))
    plot_u="1:$icol"
    ((RELATIVE)) && plot_u="1:(\$$icol/\$$icol0)"
    ((iex>0)) && plot_string="$plot_string,"
    plot_string="$plot_string '$EXSTATS' u $plot_u w l lw 1.5 ti '$iex'"
  done

  # Set y label.
  xlabel="MC step"
  ylabel="L$LNORM norm of weights"
  ((RELATIVE)) && ylabel="Relative $ylabel"

  # Dump gnuplot file.
  echo "set xlabel '$xlabel'" >> "$PLOT"
  echo "set ylabel '$ylabel'" >> "$PLOT"
  echo "$plot_string" >> "$PLOT"
  echo "pause -1" >> "$PLOT"

  # Report.
  echo "$PLOT generated.$EL"

}


main() {
  #----------------#
  # Main function. #
  #----------------#

  general_setup
  read_command_line "$@"
  parse_files
  case "$TYPE" in
  bar-video|bar-frames)
    extract_frame_data
    gen_bar_frames
    [ "$TYPE" = bar-video ] && gen_bar_video ;;
  line-plot)
    gen_plot ;;
  esac

}


########################## SCRIPT BODY ##########################


main "$@"
