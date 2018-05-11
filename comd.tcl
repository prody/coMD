# Copyright (c) 2010, University of Pittsburgh
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

package provide comd 1.0

package require solvate
package require autoionize
package require psfgen
package require autopsf
package require pbctools
package require exectool


set COMD_PATH $env(COMD_PATH)
set PACKAGE_PATH "$COMD_PATH"
set PACKAGEPATH "$COMD_PATH"

variable platform $tcl_platform(platform)
switch $platform {
  unix {
    set TMPDIR "/tmp" ;  # or even $::env(TMPDIR), at times.
  } macintosh {
    set TMPDIR $::env(TRASH_FOLDER)  ;# a better place?
  } default {
    set TMPDIR [pwd]
    catch {set TMPDIR $::env(TMP)}
    catch {set TMPDIR $::env(TEMP)}
  }
}

namespace eval ::comd:: {
  namespace export comd

  variable version 1.0

  variable w

  # Variables for system setup
  variable molid -1
  # input files
  variable walker1_pdb 
  variable walker2_pdb 
  variable walker1_chid 
  variable walker2_chid 
  # Ionization parameters
  variable topo_file [list]
  variable solvent_padding_x
  variable solvent_padding_y
  variable solvent_padding_z
  # Minimization parameters
  variable para_file [list]
  variable temperature 
  variable min_length 
  # ANM-MC-Metropolis parameters
  variable anm_cutoff 
  variable dev_mag 
  variable accept_para
  variable max_steps
  variable step_cutoff
  # TMD options
  variable spring_k 
  variable tmd_len 
  # Simulation options
  variable comd_cycle 
  variable num_cores
  variable gpus_selected
  variable gpus_selection1
  variable gpus_selection2
  variable gpus_present
  variable python_path ""
  variable NAMD_PATH ""
  # output options
  variable outputdir 
  variable output_prefix
  variable from_commandline 0
  variable run_now
  variable start_dir
  
  # Logvew window counter
  variable logcount 0
  variable lognames [list]
  variable titles [list "Prepare System"]
  variable interfaces [list "prepare"]
  variable which_mode [lindex $titles 0]
}

proc wsplit {string sep} {
    set first [string first $sep $string]
    if {$first == -1} {
        return [list $string]
    } else {
        set l [string length $sep]
        set left [string range $string 0 [expr {$first-1}]]
        set right [string range $string [expr {$first+$l}] end]
        return [concat [list $left] [wsplit $right $sep]]
    }
}

proc comd::Logview {log_file_name} {
  variable logcount
  variable lognames
  set logindex [lsearch $lognames $log_file_name]
  set log .somenonsense
  if {$logindex > -1} {
    set windowname "log$logindex"
    set log .$windowname
  }
  if {[winfo exists $log] == 0} {
    if {$logindex > -1} {
      lset lognames $logindex "somenonsense"
    }
    set logindex $logcount
    lappend lognames $log_file_name
    set windowname "log$logindex"
    set log [toplevel ".$windowname"]
    wm title $log "Logfile [lindex [file split $log_file_name] end] ($log_file_name)"
    wm resizable $log 1 1
    incr logcount

    text $log.text -bg White -bd 2 \
      -yscrollcommand ".$windowname.vscr set"
    scrollbar $log.vscr -command ".$windowname.text yview"
    pack $log.text -side left -fill both -expand 1
    pack $log.vscr -side right -fill y
  }

  $log.text configure -state normal
  #set count 0
  #set tabwidth 0
  #foreach family [lsort -dictionary [font families]] {
  #    $log.text tag configure f[incr count] -font [list $family 10]
  #    $log.text insert end ${family}:\t {} \
  #            "This is a simple sampler\n" f$count
  #    set w [font measure [$log.text cget -font] ${family}:]
  #    if {$w+5 > $tabwidth} {
  #        set tabwidth [expr {$w+5}]
  #        $log.text configure -tabs $tabwidth
  #    }
  #}
  $log.text delete 1.0 end
  set logfile [open $log_file_name "r"]
  set line ""
  while {[gets $logfile line] != -1} {
    $log.text insert end "$line\n"
  }
  close $logfile
  $log.text yview moveto 1
  $log.text configure -state disabled
}

proc ::comd::comdgui {} {
  variable w

  global env


  # If already initialized, just turn on
  if [winfo exists .comdgui] {
    wm deiconify .comdgui
    raise .comdgui
    return
  }

  # Initialize window
  set w [toplevel .comdgui]
  wm title $w "COllective Molecular Dynamics v$::comd::version"
  wm resizable $w 0 0

  # Set main frame
  set mf [frame $w.main_frame]

  # VISUALIZE results
  
  # Prepare System and Simulation Files
  set mfa [frame $mf.prepare]
  # Select input files
  set mfaif [labelframe $mfa.input_files -text "Protein structures:" -bd 2]
  # Initial PDB
  grid [button $mfaif.ini_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The walker1 protein structure should be given in the standard PDB format."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaif.ini_label -text "Initial PDB:                      " -width 21] \
    -row 1 -column 1 -sticky w
  grid [entry $mfaif.ini_path -width 47 \
      -textvariable ::comd::walker1_pdb] \
    -row 1 -column 2 -columnspan 6 -sticky ew
  grid [button $mfaif.ini_browse -text "Browse" -width 14 -pady 1 -command {
      set tempfile [tk_getOpenFile \
                    -filetypes {{"PDB files" { .pdb .PDB }} {"All files" *}}]
      if {![string equal $tempfile ""]} {
        set ::comd::walker1_pdb $tempfile
      } }] \
    -row 1 -column 8 -columnspan 3 -sticky w
       
  # Final PDB
  grid [button $mfaif.fin_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The walker2 protein structure should be given in the standard PDB format."}] \
    -row 2 -column 0 -sticky w
  grid [label $mfaif.fin_label -text "Final PDB:                        " -width 21] \
    -row 2 -column 1 -sticky w
  grid [entry $mfaif.fin_path -width 47 \
      -textvariable ::comd::walker2_pdb] \
    -row 2 -column 2 -columnspan 6 -sticky ew
  grid [button $mfaif.fin_browse -text "Browse" -width 14 -pady 1 -command {
        set tempfile [tk_getOpenFile \
          -filetypes {{"PDB files" { .pdb .PDB }} {"All files" *}}]
        if {![string equal $tempfile ""]} {
          set ::comd::walker2_pdb $tempfile
        } }] \
   -row 2 -column 8 -columnspan 3 -sticky w

  grid [button $mfaif.inich_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The chain ID for the walker1 and walker2 structure which should be in the previously imported PDB file."}] \
    -row 3 -column 0 -sticky w
  grid [label $mfaif.inich_label -text "Initial PDB chain ID:         " -width 21] \
    -row 3 -column 1 -sticky w
  grid [entry $mfaif.inich_entry -width 17 \
    -textvariable ::comd::walker1_chid] \
    -row 3 -column 2 -columnspan 3 -sticky ew  

  grid [label $mfaif.separatpr_label -width 6] \
    -row 3 -column 5 -sticky w

  grid [button $mfaif.finch_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The chain ID for the walker1 structure which should be in the previously imported PDB file."}] \
    -row 3 -column 6 -sticky w
  grid [label $mfaif.finch_label -text "Final PDB chain ID:          " -width 21] \
    -row 3 -column 7 -sticky w
  grid [entry $mfaif.finch_entry -width 17 \
    -textvariable ::comd::walker2_chid] \
    -row 3 -column 8 -columnspan 3 -sticky ew
  
    
  pack $mfaif -side top -ipadx 0 -ipady 5 -fill x -expand 1

  # Enter ionization options
  set mfaio [labelframe $mfa.ionize_options -text "Ionization parameters" -bd 2]

  #Solvation box padding and counter ions
  grid [button $mfaio.padding_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "This is the half of the walker1 distance between the protein \
and its imaginary copies under periodic boundary conditions. For systems with \
probes, the resulting padding distance will be slightly larger, due to \
constraint of preserving the ratio of 20 water molecules per probe molecule."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfaio.padding_label -text "Box padding (A):             " -width 21] \
    -row 0 -column 1 -sticky w
  grid [entry $mfaio.padding_entry_x -width 5 \
    -textvariable ::comd::solvent_padding_x] \
    -row 0 -column 2 -sticky ew
  grid [entry $mfaio.padding_entry_y -width 5 \
    -textvariable ::comd::solvent_padding_y] \
    -row 0 -column 3 -sticky ew
  grid [entry $mfaio.padding_entry_z -width 5 \
    -textvariable ::comd::solvent_padding_z] \
    -row 0 -column 4 -sticky ew

  grid [label $mfaio.separatpr_label -width 6] \
    -row 1 -column 5 -sticky w

  pack $mfaio -side top -ipadx 0 -ipady 5 -fill x -expand 1

  #Topology files
  grid [button $mfaio.topo_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "Multiple topology files can be specified that will be given to psfgen package of vmd. \
Therefore, from given static structure this will create bonds, angles and various structural elements \
based on topology parameters provided. Suggested file extension is .top but others will be accepted."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaio.topo_label -text "Topology files:                  " -width 21] \
    -row 1 -column 1 -sticky w
  grid [frame $mfaio.topo_frame] \
    -row 1 -rowspan 6 -column 2 -columnspan 6 -sticky w
  scrollbar $mfaio.topo_frame.scroll -command "$mfaio.topo_frame.list yview"
  listbox $mfaio.topo_frame.list -activestyle dotbox \
    -yscroll "$mfaio.topo_frame.scroll set" \
    -width 45 -height 6 -setgrid 1 -selectmode browse \
    -listvariable ::comd::topo_file
  frame $mfaio.topo_frame.buttons
  pack $mfaio.topo_frame.list $mfaio.topo_frame.scroll -side left -fill y -expand 1

  grid [button $mfaio.topo_add -text "Add" -width 14 -pady 1 \
        -command [namespace code {
        set tempfiles [tk_getOpenFile -multiple 1 \
          -filetypes { {{Topology files} {.top .TOP .rtf .RTF .str .STR}} {{All files} {*}} }]
        if {$tempfiles!=""} {
          foreach tempfile $tempfiles {
            if {[lsearch $::comd::topo_file $tempfile] > -1} {
              tk_messageBox -type ok -title "WARNING" \
               -message "$tempfile has already been added to the list."
            } else {
              lappend ::comd::topo_file $tempfile
            }
          }
        }
      }]] \
    -row 1 -column 8 -columnspan 3 -sticky w
  grid [button $mfaio.topo_delete -text "Remove"  -width 14 -pady 1 \
      -command [namespace code {
      foreach i [.comdgui.main_frame.prepare.ionize_options.topo_frame.list curselection] {
        .comdgui.main_frame.prepare.ionize_options.topo_frame.list delete $i
      } }]] \
    -row 2 -column 8 -columnspan 3 -sticky w

  # Enter minimization options
  set mfamo [labelframe $mfa.minimize_options -text "Minimization parameters" -bd 2]
  
  #Temperature and minimization length parameters
  grid [button $mfamo.temperature_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The temperature for molecular dynamics simulation needs to be entered. The units are in Kelvin."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfamo.temperature_label -text "Temperature (K):              " -width 21] \
    -row 0 -column 1 -sticky w
  grid [entry $mfamo.temperature_entry -width 17 \
    -textvariable ::comd::temperature] \
    -row 0 -column 2 -columnspan 3 -sticky ew

  grid [label $mfamo.separatpr_label -width 6] \
    -row 0 -column 5 -sticky w

  grid [button $mfamo.length_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "After running targeted molecular dynamics simulations the walker2 structure needs to be equilibrated into a stable state. Longer minimization creates 
        more stable structures which will guarantee users have a stable targeted molecular dynamics simulation. The units are in ps. "}] \
    -row 0 -column 6 -sticky w
  grid [label $mfamo.length_label -text "Minimization length (ps):   " -width 21] \
    -row 0 -column 7 -sticky w
  grid [entry $mfamo.length_entry -width 17 \
    -textvariable ::comd::min_length] \
    -row 0 -column 8 -columnspan 3 -sticky ew
  
  #Parameter files
  grid [button $mfamo.para_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "Multiple parameter files can be specified for the force field.
        The file should be provided in par or prm format and include necessary parameters required for NAMD."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfamo.para_label -text "Parameter files:                " -width 21] \
    -row 1 -column 1 -sticky w
  grid [frame $mfamo.para_frame] \
    -row 1 -rowspan 6 -column 2 -columnspan 6 -sticky w
  scrollbar $mfamo.para_frame.scroll -command "$mfamo.para_frame.list yview"
  listbox $mfamo.para_frame.list -activestyle dotbox \
    -yscroll "$mfamo.para_frame.scroll set" \
    -width 45 -height 6 -setgrid 1 -selectmode browse \
    -listvariable ::comd::para_file
  frame $mfamo.para_frame.buttons
  pack $mfamo.para_frame.list $mfamo.para_frame.scroll \
    -side left -fill y -expand 1

  grid [button $mfamo.para_add -text "Add" -width 14 -pady 1 \
        -command [namespace code {
        set tempfiles [tk_getOpenFile -multiple 1 \
          -filetypes { {{Parameter files} {.par .PAR .prm .PRM .str .STR}} {{All files} {*}} }]
        if {$tempfiles!=""} {
          foreach tempfile $tempfiles {
            if {[lsearch $::comd::para_file $tempfile] > -1} {
              tk_messageBox -type ok -title "WARNING" \
                -message "$tempfile has already been added to the list."
            } else {
              lappend ::comd::para_file $tempfile
            }
          }
        }
      }]] \
    -row 1 -column 8 -columnspan 3 -sticky w
  grid [button $mfamo.para_delete -text "Remove"  -width 14 -pady 1 \
      -command [namespace code {
      foreach i [.comdgui.main_frame.prepare.minimize_options.para_frame.list curselection] {
        .comdgui.main_frame.prepare.minimize_options.para_frame.list delete $i
      } }]] \
    -row 2 -column 8 -columnspan 3 -sticky w

  pack $mfamo -side top -ipadx 0 -ipady 5 -fill x -expand 1

  set mfamc [labelframe $mfa.anmmc_options -text "ANM-MC-Metropolis options:" -bd 2]
  
  #grid [button $mfamc.anmcut_help -text "?" -width 1 -padx 0 -pady 0 -command {
  #    tk_messageBox -type ok -title "HELP" \
  #      -message "In ANM calculations, the cutoff parameter is the maximum distance that two residues are in contact. The units are A. "}] \
  #  -row 0 -column 0 -sticky w
  #grid [label $mfamc.anmc_label -text "ANM cutoff (A):              " -width 21] \
  #  -row 0 -column 1 -sticky w
  #grid [entry $mfamc.anmc_field -width 17 \
  #  -textvariable ::comd::anm_cutoff] \
  #  -row 0 -column 2 -columnspan 3 -sticky w

  grid [button $mfamc.step_cutoffut_help -text "?" -padx 0 -pady 0 -command {
    tk_messageBox -type ok -title "HELP" \
      -message "To keep structure intact and to avoid having unrealistic \
        and very different structures in ANM-MC step, an rmsd threshold is used. Suggested value is 4 A."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfamc.step_cutoff_label -text "Step cutoff (A): "] \
    -row 0 -column 1 -sticky w
  grid [entry $mfamc.step_cutoff_field -width 17 \
    -textvariable ::comd::step_cutoff] \
    -row 0 -column 2 -columnspan 3 -sticky w

  grid [label $mfamc.separatpr1_label -width 6] \
    -row 0 -column 5 -sticky w

  grid [button $mfamc.dev_mag_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The scaling factor used when disturbing the protein structure in ANM-MC steps. Default and suggested value is 0.1 A."}] \
    -row 0 -column 6 -sticky w
  grid [label $mfamc.dev_mag_label -text "Deviation (A):                   " -width 21] \
    -row 0 -column 7 -sticky w
  grid [entry $mfamc.dev_mag_field -width 17 \
      -textvariable ::comd::dev_mag] \
    -row 0 -column 8 -columnspan 3 -sticky w

  grid [button $mfamc.accept_para_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The starting value for the acceptance parameter in ANM-MC steps. Default and suggested value is 0.1."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfamc.accept_para_label -text "Acceptance parameter:     " -width 21] \
    -row 1 -column 1 -sticky w
  grid [entry $mfamc.accept_para_field -width 17 \
      -textvariable ::comd::accept_para] \
    -row 1 -column 2 -columnspan 3 -sticky w

  grid [label $mfamc.separatpr2_label -width 6] \
    -row 1 -column 5 -sticky w

  grid [button $mfamc.max_steps_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The maximal number of steps in ANM-MC step. Default and suggested value is 1000000."}] \
    -row 1 -column 6 -sticky w
  grid [label $mfamc.max_steps_label -text "Max no of ANM steps:      " -width 21] \
    -row 1 -column 7 -sticky w
  grid [entry $mfamc.max_steps_field -width 17 \
      -textvariable ::comd::max_steps] \
    -row 1 -column 8 -columnspan 3 -sticky w

  pack $mfamc -side top -ipadx 0 -ipady 5 -fill x -expand 1

  set mfatm [labelframe $mfa.tmd_options -text "TMD options:" -bd 2]

  grid [button $mfatm.spring_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "In targeted molecular dynamics simulation, the target potential is harmonic and \
the spring constant term shows the force applied to a given structure to reach the target structure."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfatm.spring_label -text "Spring constant:              " -width 21] \
    -row 0 -column 1 -sticky w
  grid [entry $mfatm.spring_field -width 17 \
      -textvariable ::comd::spring_k] \
    -row 0 -column 2 -columnspan 3 -sticky w

  grid [label $mfatm.separatpr2_label -width 6] \
    -row 0 -column 5 -sticky w

  grid [button $mfatm.tmd_len_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The length of targeted molecular dynamics simulations in the units of ps. The length of collective molecular dynamics will change based on the structure and for a structure with 200 residues suggested length is in the order of hundreds."}] \
    -row 0 -column 6 -sticky w
  grid [label $mfatm.tmd_len_label -text "TMD length (ps):              " -width 21] \
    -row 0 -column 7 -sticky w
  grid [entry $mfatm.tmd_len_field -width 17 \
      -textvariable ::comd::tmd_len] \
    -row 0 -column 8 -columnspan 3 -sticky w
  
  pack $mfatm -side top -ipadx 0 -ipady 5 -fill x -expand 1

  ########################################3
  set mfaso [labelframe $mfa.simulation_options -text "Simulation options:" -bd 2]

  grid [button $mfaso.cmd_cyc_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "Each CoMD cycle consists of minimization, ANM-MC-Metropolis disturbances and targeted molecular dynamics. Please choose the maximum number of cycles performed. Fewer cycles may be run if the starting and walker2 structures for a given cycle are very close."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfaso.cmd_cyc_label -text "No of coMD cycles:         " -width 21] \
    -row 0 -column 1 -sticky w
  grid [entry $mfaso.cmd_cyc_field -width 17 \
      -textvariable ::comd::comd_cycle] \
    -row 0 -column 2 -columnspan 3 -sticky w

  grid [label $mfaso.separatpr1_label -width 6] \
    -row 0 -column 5 -sticky w

  grid [button $mfaso.run_now_help -text "?" -width 1 -padx 0 -pady 0 -command {
    tk_messageBox -type ok -title "HELP" \
      -message "If this is checked the simulations will run as soon as the system is prepared"}] \
    -row 0 -column 6 -sticky w
  grid [label $mfaso.run_now_label -text "Run now:                       " -width 21] \
    -row 0 -column 7 -sticky w
  grid [label $mfaso.separatpr2_label -width 13] \
    -row 1 -column 8 -columnspan 2 -sticky w
  grid [checkbutton $mfaso.run_now_check -width 1 \
      -variable ::comd::run_now] \
    -row 0 -column 10 -sticky e

  grid [button $mfaso.gpu_id_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The identifiers for the GPUs that will run your TMD simulation separated by commas. NAMD can use one GPU per thread and multiple threads can share GPUs."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaso.gpu_id_label -text "GPU IDs:                        " -width 21] \
    -row 1 -column 1 -sticky w
  grid [entry $mfaso.gpu_id_field -width 17 \
      -textvariable ::comd::gpus_selected] \
    -row 1 -column 2 -columnspan 3 -sticky w

  grid [label $mfaso.separatpr3_label -width 6] \
    -row 1 -column 5 -sticky w

  grid [button $mfaso.num_cores_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The number of physical cores in the cluster or PC that will run your TMD simulation. NAMD is running parallel on CPUs."}] \
    -row 1 -column 6 -sticky w
  grid [label $mfaso.num_cores_label -text "No of physical cores:         " -width 21] \
    -row 1 -column 7 -sticky w
  grid [entry $mfaso.num_cores -width 17 \
      -textvariable ::comd::num_cores] \
    -row 1 -column 8 -columnspan 3 -sticky ew

  
  pack $mfaso -side top -ipadx 0 -ipady 5 -fill x -expand 1



  ########################################4


  set mfaoo [labelframe $mfa.output_options -text "Output options:" -bd 2]

  grid [button $mfaoo.outdir_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "Output folder, default is current working directory."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfaoo.outdir_label -text "Output folder:                  " -width 21] \
    -row 0 -column 1 -sticky w
  grid [entry $mfaoo.outdir_path -width 47 -textvariable ::comd::outputdir] \
    -row 0 -column 2 -columnspan 6 -sticky ew
  grid [button $mfaoo.dcd_browse -text "Browse" -width 14 -pady 1 -command {
      set tempfile [tk_chooseDirectory]
      if {![string equal $tempfile ""]} {
        set ::comd::outputdir $tempfile
      }}] \
    -row 0 -column 8 -columnspan 3 -sticky w

  grid [button $mfaoo.prefix_help -text "?" -width 1 -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "All output files and folders will start with this prefix.\
A unique and descriptive prefix choice may allow running multiple simulations in the same folder."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaoo.prefix_label -text "Output prefix:                  " -width 21] \
    -row 1 -column 1 -sticky w
  grid [entry $mfaoo.prefix_path -width 17 \
      -textvariable ::comd::output_prefix] \
    -row 1 -column 2 -columnspan 3 -sticky w

  grid [label $mfaoo.separator_label -width 6] \
    -row 1 -column 5 -sticky w
  
  pack $mfaoo -side top -ipadx 0 -ipady 5 -fill x -expand 1

  # Prepare System
  button $mfa.button -text "Prepare System" -command ::comd::Prepare_system -bd 3
  pack $mfa.button

  pack $mfa -side top -padx 0 -pady 0 -fill x -expand 1
  pack $mf -side top -padx 0 -pady 0 -fill x -expand 1

  return $w
}

proc ::comd::Prepare_system {} {

  # WHAT IS NEW?
  # 2.1 - Bug fixes, and file checks
  # 2.0 - Improved system setup provides lesser number of solvent atoms
  # 2.0 - Cleans up intermediate files
  # 2.0 - Outputs a log file for troubleshooting, and further intstructions
  # 2.0 - NAMD configuration files are prepared for a single or multiple
  #       simulations

  # HOW THE CODE WORKS
  # The code will
  # (1)   solvate the protein, or everything in the PDB/PSF files that you provide
  # (2)   add counter ions to neutralize the system
  # (3)   write output files for each simulation

  #       Setups of multiple simulations differ only at random number seeds.
  #       This will be sufficient to result in a different trajectory.
  variable w
  variable pro

  if {$::comd::outputdir != ""} {
      if {![file isdirectory $::comd::outputdir]} {
        if {[catch {file mkdir $::comd::outputdir}]} {
          if {[info exists ::comd::from_commandline]} {
            error "Could not make output folder: $::comd::outputdir"
          } else {
            tk_messageBox -type ok -title "ERROR" \
            -message "Could not make output folder: $::comd::outputdir"
          }
          
          return

        }
      }
  }

  if {$::comd::walker1_pdb == "" || $::comd::walker2_pdb == ""} {
    tk_messageBox -type ok -title "ERROR" \
      -message "Both PDB files must be specified."
    return
  }

  if {$::comd::solvent_padding_z < 4} {
    tk_messageBox -type ok -title "ERROR" \
      -message "Solvent box padding parameter must be larger than 4 A at least in the z direction."
    return
  }

  if {[string length [string trim $::comd::output_prefix]] == 0} {
    tk_messageBox -type ok -title "ERROR" \
      -message "Please enter a descriptive name (prefix) for the system."
    return
  }

  if {$::comd::temperature == ""} {
    tk_messageBox -type ok -title "ERROR" \
      -message "The temperature must be specified."
    return
  }

  global env
  global COMD_PATH

  set log_file [open [file join "$::comd::outputdir" "$::comd::output_prefix.log"] w]
  puts $log_file "---==## [clock format [clock seconds]] #==---"
  puts $log_file "Version: $::comd::version"
  puts $log_file "Info: Logging started for setup of $::comd::output_prefix."
  puts $log_file "Solvation: Box padding $::comd::solvent_padding_x A in x, $::comd::solvent_padding_y A in y, $::comd::solvent_padding_z A in z."

  
  ####### SOLVATION AND IONIZATION OF WALKER 1 PROTEIN STRUCTURE #######
  resetpsf
  mol delete all
  mol new $::comd::walker1_pdb
  if {[info exists ::comd::walker1_chid]} {
    set pro [atomselect top "not altloc B and not hydrogen and chain $::comd::walker1_chid"]
  } else {
    set pro [atomselect top "not altloc B and not hydrogen and protein and not resname UNK"]
  }
  
  $pro writepdb init.pdb
  mol delete all
  mol new init.pdb

  if {$::comd::topo_file == [list]} {
    lappend ::comd::topo_file "$COMD_PATH/top_all36_prot.rtf"
    lappend ::comd::topo_file "$COMD_PATH/toppar_water_ions.str"
  }
  autopsf -mol top -top $::comd::topo_file -prefix pro
  solvate pro_formatted_autopsf.psf pro_formatted_autopsf.pdb \
    -x $::comd::solvent_padding_x -y $::comd::solvent_padding_y -z $::comd::solvent_padding_z \
    +x $::comd::solvent_padding_x +y $::comd::solvent_padding_y +z $::comd::solvent_padding_z -o pro_wb

  set totalcharge 0
  foreach charge [[atomselect top "all"] get charge] {
    set totalcharge [expr $totalcharge + $charge]
  }
  # number of CL and NA atoms are determined
  set nna 0
  set ncl 0
  puts $log_file "Ionization: Initial PDB System has a total charge of $totalcharge electrons."
  if {$totalcharge > 0} {
      set ncl [expr round($totalcharge)]
      puts $log_file "Ionization: $ncl chloride ions will be added to walker1."
      autoionize -psf pro_wb.psf -pdb pro_wb.pdb \
      -o [file join $::comd::outputdir "walker1_ionized"] -from 5 -between 5 -nna $nna -ncl $ncl -seg ION
      puts $log_file "Ionization: Initial PDB System is ionized to become neutral."
  } elseif {$totalcharge < 0} {
      set nna [expr -1 * round($totalcharge)]
      puts $log_file "Ionization: $nna sodium ions will be added to walker1."
      autoionize -psf pro_wb.psf -pdb pro_wb.pdb \
      -o [file join $::comd::outputdir "walker1_ionized"] -from 5 -between 5 -nna $nna -ncl $ncl -seg ION
      puts $log_file "Ionization: Initial PDB System is ionized to become neutral."
  }
  
  mol new init.pdb
  set everyone [atomselect top "all"]
  set walker1_lims [measure minmax $everyone]
  set walker1_cent [measure center $everyone]
  set ixmin [lindex $walker1_lims 0 0]
  set ixmax [lindex $walker1_lims 1 0]
  set iymin [lindex $walker1_lims 0 1]
  set iymax [lindex $walker1_lims 1 1]
  set izmin [lindex $walker1_lims 0 2]
  set izmax [lindex $walker1_lims 1 2]
  set ixcen [lindex $walker1_cent 0]
  set iycen [lindex $walker1_cent 1]
  set izcen [lindex $walker1_cent 2]
  set ixlen [expr {$ixmax-$ixmin+16}]
  set iylen [expr {$iymax-$iymin+16}]
  set izlen [expr {$izmax-$izmin+16}]

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    ####### SOLVATION AND IONIZATION OF WALKER 2 PROTEIN STRUCTURE #######
    resetpsf
    mol delete all
    mol new $::comd::walker2_pdb
    if {[info exists ::comd::walker2_chid]} {
      set pro [atomselect top "not altloc B and not hydrogen and chain $::comd::walker2_chid"]
    } else {
      set pro [atomselect top "not altloc B and not hydrogen and protein and not resname UNK"]
    }
    $pro writepdb fino.pdb
    mol delete all
    mol new fino.pdb
    autopsf -mol top -top $::comd::topo_file -prefix pro
    solvate pro_formatted_autopsf.psf pro_formatted_autopsf.pdb \
      -x $::comd::solvent_padding_x -y $::comd::solvent_padding_y -z $::comd::solvent_padding_z \
      +x $::comd::solvent_padding_x +y $::comd::solvent_padding_y +z $::comd::solvent_padding_z -o pro_wb

    set totalcharge 0
    foreach charge [[atomselect top "all"] get charge] {
      set totalcharge [expr $totalcharge + $charge]
    }
    # number of CL and NA atoms are determined
    set nna 0
    set ncl 0
    puts $log_file "Ionization: Final PDB System has a total charge of $totalcharge electrons."
    if {$totalcharge > 0} {
      set ncl [expr round($totalcharge)]
      puts $log_file "Ionization: $ncl chloride ions will be added to walker2."
      autoionize -psf pro_wb.psf -pdb pro_wb.pdb \
      -o [file join $::comd::outputdir "walker2_ionized"] -from 5 -between 5 -nna $nna -ncl $ncl -seg ION
      puts $log_file "Ionization: Final PDB System is ionized to become neutral."
    } elseif {$totalcharge < 0} {
      set nna [expr -1 * round($totalcharge)]
      puts $log_file "Ionization: $nna sodium ions will be added to walker2."
      autoionize -psf pro_wb.psf -pdb pro_wb.pdb \
      -o [file join $::comd::outputdir "walker2_ionized"] -from 5 -between 5 -nna $nna -ncl $ncl -seg ION
      puts $log_file "Ionization: Final PDB System is ionized to become neutral."
    }

    mol new fino.pdb
    set everyone [atomselect top "all"]
    set walker2_lims [measure minmax $everyone]
    set walker2_cent [measure center $everyone]
    set fxmin [lindex $walker2_lims 0 0]
    set fxmax [lindex $walker2_lims 1 0]
    set fymin [lindex $walker2_lims 0 1]
    set fymax [lindex $walker2_lims 1 1]
    set fzmin [lindex $walker2_lims 0 2]
    set fzmax [lindex $walker2_lims 1 2]
    set fxcen [lindex $walker2_cent 0]
    set fycen [lindex $walker2_cent 1]
    set fzcen [lindex $walker2_cent 2]
    set fxlen [expr {$fxmax-$fxmin+16.0}]
    set fylen [expr {$fymax-$fymin+16.0}]
    set fzlen [expr {$fzmax-$fzmin+16.0}]
  }
  
  ####### INITIAL MINIMIZATION OF STARTING PROTEIN STRUCTURES #######
  puts $log_file "Simulation: NAMD configuration files for minimization written in ${::comd::output_prefix}_walker1_min"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $log_file "and ${::comd::output_prefix}_walker2_min."
  }
  close $log_file

  if {$::comd::para_file == [list]} {
    lappend ::comd::para_file "$COMD_PATH/par_all36_prot.prm" 
    lappend ::comd::para_file "$COMD_PATH/par_all36m_prot.prm" 
    lappend ::comd::para_file "$COMD_PATH/toppar_water_ions.str"
  }
  set tcl_file_name [file join "$::comd::outputdir" "$::comd::output_prefix.tcl"]
  set tcl_file [open $tcl_file_name w] 
  puts $tcl_file "#This tcl file will run full collective molecular dynamics simulation with given parameters."
  puts $tcl_file "cd $::comd::outputdir"
  puts $tcl_file "set sh_filename \"${::comd::output_prefix}_min0.sh\""
  puts $tcl_file "set sh_file \[open \$sh_filename w\]"

  puts $tcl_file "package require exectool"
  if {$::comd::NAMD_PATH == ""} {
    puts $tcl_file "set namd2path \[::ExecTool::find \"namd2\"\]"
  } else {
    puts $tcl_file "set namd2path ${::comd::NAMD_PATH}"
  }

  if {$::comd::python_path == ""} {
    puts $tcl_file "set python_path \[::ExecTool::find \"python\"\]"	
  } else {
    puts $tcl_file "set python_path $::comd::python_path\/python" 
  }


  puts $tcl_file "puts \$sh_file \"\\\#\\\!\\\/bin\\\/bash\""

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    set ::comd::num_cores [expr {$::comd::num_cores / 2}]
  }

  if {$::comd::gpus_present} {
    set processes_per_run [expr {[llength [wsplit $::comd::gpus_selection1 ","]] + 1}]

    if {[info exists ::comd::num_cores]} {
      set remainder [expr {$::comd::num_cores % $processes_per_run}]
      set processes_per_run [expr {$::comd::num_cores - $remainder - $processes_per_run}]
    }
  } else {
    set processes_per_run [expr {$::comd::num_cores - 1}] 
  }

  puts $tcl_file "puts \$sh_file \"NAMD=\\\"\$namd2path \+idlepoll \\\"\""

  # Walker 1 minimization
  puts $tcl_file "file mkdir \"${::comd::output_prefix}_walker1_min\""
  puts $tcl_file "set namd_file \[open \[file join \"${::comd::output_prefix}_walker1_min\" \"min.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ..\/walker1_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ..\/walker1_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $::comd::temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname walker1_minimized0\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""

  foreach tempfile $::comd::para_file {
    puts $tcl_file "puts \$namd_file \"parameters $tempfile\""
  }

  puts $tcl_file "puts \$namd_file \"temperature $::comd::temperature\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector1 ${ixlen},0,0\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector2 0,${iylen},0\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector3 0,0,${izlen}\""
  puts $tcl_file "puts \$namd_file \"cellOrigin ${ixcen},${iycen},${izcen}\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies $::comd::min_length\""
  puts $tcl_file "puts \$namd_file \"outputPressure $::comd::min_length\""
  puts $tcl_file "puts \$namd_file \"restartfreq $::comd::min_length\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr $::comd::min_length*5]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr $::comd::min_length*5]\""
  puts $tcl_file "puts \$namd_file \"minimize [expr $::comd::min_length*5]\""
  puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${::comd::output_prefix}_walker1_min\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD \+devices $::comd::gpus_selection1 \+ppn $processes_per_run min.conf > min0.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\"" 

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    # Walker 2 minimization
    puts $tcl_file "file mkdir \"${::comd::output_prefix}_walker2_min\""
    puts $tcl_file "set namd_file \[open \[file join \"${::comd::output_prefix}_walker2_min\" \"min.conf\"\] w\]"
    puts $tcl_file "puts \$namd_file \"coordinates     ..\/walker2_ionized.pdb\""
    puts $tcl_file "puts \$namd_file \"structure       ..\/walker2_ionized.psf\""
    puts $tcl_file "puts \$namd_file \"set temperature $::comd::temperature\""
    puts $tcl_file "puts \$namd_file \"set outputname walker2_minimized0\""
    puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
    puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""

    foreach tempfile $::comd::para_file {
      puts $tcl_file "puts \$namd_file \"parameters $tempfile\""
    }

    puts $tcl_file "puts \$namd_file \"temperature $::comd::temperature\""
    puts $tcl_file "puts \$namd_file \"cellBasisVector1 ${fxlen},0,0\""
    puts $tcl_file "puts \$namd_file \"cellBasisVector2 0,${fylen},0\""
    puts $tcl_file "puts \$namd_file \"cellBasisVector3 0,0,${fzlen}\""
    puts $tcl_file "puts \$namd_file \"cellOrigin ${fxcen},${fycen},${fzcen}\""
    puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
    puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
    puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
    puts $tcl_file "puts \$namd_file \"timestep        1.0\""
    puts $tcl_file "puts \$namd_file \"switching       on\""
    puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
    puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
    puts $tcl_file "puts \$namd_file \"rigidBonds none\""
    puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
    puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
    puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
    puts $tcl_file "puts \$namd_file \"langevin on\""
    puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
    puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
    puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
    puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
    puts $tcl_file "puts \$namd_file \"outputEnergies $::comd::min_length\""
    puts $tcl_file "puts \$namd_file \"outputPressure $::comd::min_length\""
    puts $tcl_file "puts \$namd_file \"restartfreq $::comd::min_length\""
    puts $tcl_file "puts \$namd_file \"dcdfreq [expr $::comd::min_length*5]\""
    puts $tcl_file "puts \$namd_file \"xstfreq [expr $::comd::min_length*5]\""
    puts $tcl_file "puts \$namd_file \"minimize [expr $::comd::min_length*5]\""
    puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
    puts $tcl_file "close \$namd_file"
    puts $tcl_file "puts \$sh_file \"cd ${::comd::output_prefix}_walker2_min\""
    puts $tcl_file "puts \$sh_file \"\\\$NAMD \+devices $::comd::gpus_selection2 \+ppn $processes_per_run min.conf > min0.log \&\""
    puts $tcl_file "puts \$sh_file \"cd ..\"" 
  }

  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "puts \"Now running minimization 0\""
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"

  puts $tcl_file "if {\$status} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: \$output\""
  puts $tcl_file "exit"
  puts $tcl_file "}"

  puts $tcl_file "puts \"Finished minimization 0\""

  puts $tcl_file "set status \[catch \{exec cp ${::comd::output_prefix}_walker1_min\/walker1_minimized0.dcd initr.dcd\} output\]"
  puts $tcl_file "if {\$status} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: \$output\""
  puts $tcl_file "exit"
  puts $tcl_file "}"
  
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "set status \[catch \{exec cp ${::comd::output_prefix}_walker2_min\/walker2_minimized0.dcd fintr.dcd\} output\]" 
    puts $tcl_file "if {\$status} {"
    puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
    puts $tcl_file "puts \$err_file \"ERROR: \$output\""
    puts $tcl_file "exit"
    puts $tcl_file "}"
  }

  puts $tcl_file "puts \"Finished copying step after minimization 0\""

  # Check if any files are missing and if so raise an error
  puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker1_min/walker1_minimized0.coor r} fid\]} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: The original minimization of structure 1 failed. See ${::comd::output_prefix}_walker1_min/min0.log for details.\""
  puts $tcl_file "exit"
  puts $tcl_file "}"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker2_min/walker2_minimized0.coor r} fid\]} {"
    puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
    puts $tcl_file "puts \$err_file \"ERROR: The original minimization of structure 2 failed. See ${::comd::output_prefix}_walker2_min/min0.log for details.\""
    puts $tcl_file "exit"
    puts $tcl_file "}"
  }

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    #puts $tcl_file "package require psfgen"
    puts $tcl_file "mol delete all" 
    puts $tcl_file "mol load psf walker1_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker1_min/walker1_minimized0.coor" 
    puts $tcl_file "set sel1 \[atomselect top \"name CA\"\]" 
    puts $tcl_file "set sel1a \[atomselect top all\]"
    puts $tcl_file "mol load psf walker2_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker2_min/walker2_minimized0.coor"  
    puts $tcl_file "set sel2 \[atomselect top \"name CA\"\]" 
    puts $tcl_file "set sel2a \[atomselect top all\]"
    puts $tcl_file "set trans_mat \[measure fit \$sel2 \$sel1\]"
    puts $tcl_file "\$sel2a move \$trans_mat"
    puts $tcl_file "set rmsd \[measure rmsd \$sel2 \$sel1\]"
    puts $tcl_file "set all_rmsd(0) \$rmsd"
    puts $tcl_file "set rmsd_filename rmsd.txt"
    puts $tcl_file "set rmsd_file \[open \$rmsd_filename w\]"
    puts $tcl_file "puts \$rmsd_file \"\$rmsd\""
    puts $tcl_file "file mkdir ${::comd::output_prefix}_walker2_pro"
  }

  puts $tcl_file "file mkdir ${::comd::output_prefix}_walker1_pro"

  #loop start
  puts $tcl_file "for {set cycle 1} {\$cycle < $::comd::comd_cycle} {incr cycle} {"
  puts $tcl_file "mol delete all"

  # Check if any files are missing and if so retry the previous cycle
  puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker1_min/walker1_minimized\[expr \$\{cycle\}-1\].coor r} fid\]} {"
  puts $tcl_file "set cycle \[expr \$\{cycle\}-1\]"
  puts $tcl_file "}"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker2_min/walker2_minimized\[expr \$\{cycle\}-1\].coor r} fid\]} {"
    puts $tcl_file "set cycle \[expr \$\{cycle\}-1\]"
    puts $tcl_file "}"
  }

  # Prepare PDB files for ANM-MC (making sure they are aligned if transitioning)
  puts $tcl_file "mol load psf walker1_ionized.psf"
  puts $tcl_file "mol addfile ${::comd::output_prefix}_walker1_min/walker1_minimized\[expr \$\{cycle\}-1\].coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "set s2 \[atomselect top \"all\"\]"
    puts $tcl_file "mol load psf walker2_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker2_min/walker2_minimized\[expr \$\{cycle\}-1\].coor"
    puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
    puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
    puts $tcl_file "\$s2 move \$trans_mat"
  }
  puts $tcl_file "\$s1 writepdb starting_walker1.pdb"

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "mol delete all"
    puts $tcl_file "mol load psf walker2_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker2_min/walker2_minimized\[expr \$\{cycle\}-1\].coor"
    puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
    puts $tcl_file "set s2 \[atomselect top \"all\"\]"
    puts $tcl_file "mol load psf walker1_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker1_min/walker1_minimized\[expr \$\{cycle\}-1\].coor"
    puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
    puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
    puts $tcl_file "\$s2 move \$trans_mat"
    puts $tcl_file "\$s1 writepdb starting_walker2.pdb"
    
    puts $tcl_file "mol delete all"
    puts $tcl_file "mol load psf walker2_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker2_min/walker2_minimized\[expr \$\{cycle\}-1\].coor"
    puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  }
  puts $tcl_file "\$s1 writepdb walker1_target.pdb"

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "mol delete all"
  
    puts $tcl_file "mol load psf walker1_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker1_min/walker1_minimized\[expr \$\{cycle\}-1\].coor"
    puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
    puts $tcl_file "\$s1 writepdb walker2_target.pdb"
  }
  
  if {$::comd::anm_cutoff eq ""} {set ::comd::anm_cutoff 0}
  if {$::comd::accept_para eq ""} {set ::comd::accept_para 0}
  if {$::comd::max_steps eq ""} {set ::comd::max_steps 0}
  if {$::comd::step_cutoff eq ""} {set ::comd::step_cutoff 0}
  puts $tcl_file "set sh_file \[open \"${::comd::output_prefix}_anmmc_\$cycle.sh\" w\]"
  puts $tcl_file "set sh_filename \"${::comd::output_prefix}_anmmc_\$cycle.sh\""
  if {[info exists ::comd::num_cores]} {
    puts $tcl_file "puts \$sh_file \"export MKL_NUM_THREADS=$::comd::num_cores\""
  }
  puts $tcl_file "puts \$sh_file \"\$python_path anmmc.py starting_walker1.pdb \
    walker1_target.pdb $::comd::walker1_pdb $::comd::walker2_pdb \$cycle \$::comd::dev_mag \
    \$::comd::step_cutoff \$::comd::accept_para \$::comd::anm_cutoff \$::comd::max_steps \
    \>& cycle_\${cycle}_ini_anmmc_log.txt \&\""
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "puts \$sh_file \"\$python_path anmmc.py starting_walker2.pdb \
    walker2_target.pdb $::comd::walker1_pdb $::comd::walker2_pdb \$cycle \$::comd::dev_mag \
    \$::comd::step_cutoff \$::comd::accept_para \$::comd::anm_cutoff \$::comd::max_steps \
    \>& cycle_\${cycle}_fin_anmmc_log.txt \&\""
  }
  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "puts \"Now running ANM Monte-Carlo stepping \$cycle\""
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"
  puts $tcl_file "puts \"Finished ANM Monte-Carlo stepping \$cycle\""

  # Check if any files are missing and if so raise an error
  puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker1_min/walker1_minimized0.coor r} fid\]} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: The original minimization of structure 1 failed. Please try again with a different structure 1.\""
  puts $tcl_file "exit"
  puts $tcl_file "}"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker2_min/walker2_minimized0.coor r} fid\]} {"
    puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
    puts $tcl_file "puts \$err_file \"ERROR: The original minimization of structure 2 failed. Please try again with a different structure 2.\""
    puts $tcl_file "exit"
    puts $tcl_file "}"
  }
  
  # Prepare adjusted Calpha positions as inputs for TMD (aligning if transitioning)
  # Check for missing files and raise errors along the way.
  puts $tcl_file "mol delete all"
  puts $tcl_file "mol load psf walker1_ionized.psf"
  puts $tcl_file "mol addfile ${::comd::output_prefix}_walker1_min/walker1_minimized\[expr \$\{cycle\}-1\].coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set s2 \[atomselect top \"all\"\]"

  puts $tcl_file "mol load pdb starting_walker1.pdb"
  puts $tcl_file "if {\[catch {mol addfile cycle_\$\{cycle\}_starting_walker1_walker1_target_final_structure.dcd} \]} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: ANM-MC stepping for walker 1 (ini) in cycle \$cycle failed. See the relevant log file.\""
  puts $tcl_file "exit"
  puts $tcl_file "}"

  puts $tcl_file "mol addfile cycle_\$\{cycle\}_starting_walker1_walker1_target_final_structure.dcd"
  puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
    puts $tcl_file "\$s3 move \$trans_mat"
  }

  puts $tcl_file "\$s1 set \{x y z\} \[\$s3 get \{x y z\}\]"
  puts $tcl_file "\$s2 set occupancy 0"
  puts $tcl_file "\$s1 set occupancy 1"
  puts $tcl_file "\$s2 writepdb walker1_adjust.pdb"

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "mol delete all"
    puts $tcl_file "mol load psf walker2_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker2_min/walker2_minimized\[expr \$\{cycle\}-1\].coor"
    puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
    puts $tcl_file "set s2 \[atomselect top \"all\"\]"

    puts $tcl_file "mol load pdb starting_walker2.pdb"
    puts $tcl_file "if {\[catch {mol addfile cycle_\$\{cycle\}_starting_walker2_walker2_target_final_structure.dcd} \]} {"
    puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
    puts $tcl_file "puts \$err_file \"ERROR: ANM-MC stepping for walker 2 (ini) in cycle \$cycle failed. See the relevant log file.\""
    puts $tcl_file "exit"
    puts $tcl_file "}"
 
    puts $tcl_file "mol addfile cycle_\$\{cycle\}_starting_walker2_walker2_target_final_structure.dcd"
    puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
    puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
    puts $tcl_file "\$s3 move \$trans_mat"
    puts $tcl_file "\$s1 set \{x y z\} \[\$s3 get \{x y z\}\]"
    puts $tcl_file "\$s2 set occupancy 0"
    puts $tcl_file "\$s1 set occupancy 1"
    puts $tcl_file "\$s2 writepdb walker2_adjust.pdb"
  }

  ###### TARGETED MD ########
  puts $tcl_file "set sh_file \[open \"${::comd::output_prefix}_tmd_\$cycle.sh\" w\]"
  puts $tcl_file "set sh_filename \"${::comd::output_prefix}_tmd_\$cycle.sh\""
  puts $tcl_file "puts \$sh_file \"\\\#\\\!\\\/bin\\\/bash\""
  puts $tcl_file "puts \$sh_file \"NAMD=\\\"\$namd2path \\\"\""

  # Walker 1 TMD
  puts $tcl_file "set namd_file \[open \[file join \"${::comd::output_prefix}_walker1_pro\" \"pro.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ..\/walker1_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ..\/walker1_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $::comd::temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname walker1_process\$\{cycle\}\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""

  foreach tempfile $::comd::para_file {
    puts $tcl_file "puts \$namd_file \"parameters $tempfile\""
  }

  puts $tcl_file "puts \$namd_file \"set restartname res\""
  puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${::comd::output_prefix}_walker1_min\/walker1_minimized\[expr \$\{cycle\}-1\].coor\""
  puts $tcl_file "puts \$namd_file \"binvelocities ..\/${::comd::output_prefix}_walker1_min\/walker1_minimized\[expr \$\{cycle\}-1\].vel\""
  puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${::comd::output_prefix}_walker1_min\/walker1_minimized\[expr \$\{cycle\}-1\].xst\""
  puts $tcl_file "puts \$namd_file \"wrapWater on\""
  puts $tcl_file "puts \$namd_file \"wrapAll on\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"PME yes\""
  puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"TMD on\""
  puts $tcl_file "puts \$namd_file \"TMDk $::comd::spring_k\""
  puts $tcl_file "puts \$namd_file \"TMDOutputFreq [expr $::comd::tmd_len*1]\""
  puts $tcl_file "puts \$namd_file \"TMDFile ..\/walker1_adjust.pdb\""
  puts $tcl_file "puts \$namd_file \"TMDFirstStep 0\""
  puts $tcl_file "puts \$namd_file \"TMDLastStep [expr $::comd::tmd_len*5]\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies [expr $::comd::tmd_len*1]\""
  puts $tcl_file "puts \$namd_file \"outputPressure [expr $::comd::tmd_len*1]\""
  puts $tcl_file "puts \$namd_file \"restartfreq [expr $::comd::tmd_len*1]\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr $::comd::tmd_len*1]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr $::comd::tmd_len*1]\""
  puts $tcl_file "puts \$namd_file \"run [expr $::comd::tmd_len*5]\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${::comd::output_prefix}_walker1_pro\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD \+devices $::comd::gpus_selection1 \+ppn $processes_per_run pro.conf > pro\$\{cycle\}.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\""

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    # Walker 2 TMD
    puts $tcl_file "set namd_file \[open \[file join \"${::comd::output_prefix}_walker2_pro\" \"pro.conf\"\] w\]"
    puts $tcl_file "puts \$namd_file \"coordinates     ../walker2_ionized.pdb\""
    puts $tcl_file "puts \$namd_file \"structure       ../walker2_ionized.psf\""
    puts $tcl_file "puts \$namd_file \"set temperature $::comd::temperature\""
    puts $tcl_file "puts \$namd_file \"set outputname walker2_process\$\{cycle\}\""
    puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
    puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""

    foreach tempfile $::comd::para_file {
      puts $tcl_file "puts \$namd_file \"parameters $tempfile\""
    }

    puts $tcl_file "puts \$namd_file \"set restartname res\""
    puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${::comd::output_prefix}_walker2_min\/walker2_minimized\[expr \$\{cycle\}-1\].coor\""
    puts $tcl_file "puts \$namd_file \"binvelocities ..\/${::comd::output_prefix}_walker2_min\/walker2_minimized\[expr \$\{cycle\}-1\].vel\""
    puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${::comd::output_prefix}_walker2_min\/walker2_minimized\[expr \$\{cycle\}-1\].xst\""
    puts $tcl_file "puts \$namd_file \"wrapWater on\""
    puts $tcl_file "puts \$namd_file \"wrapAll on\""
    puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
    puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
    puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
    puts $tcl_file "puts \$namd_file \"timestep        1.0\""
    puts $tcl_file "puts \$namd_file \"switching       on\""
    puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
    puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
    puts $tcl_file "puts \$namd_file \"rigidBonds none\""
    puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
    puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
    puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
    puts $tcl_file "puts \$namd_file \"PME yes\""
    puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
    puts $tcl_file "puts \$namd_file \"langevin on\""
    puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
    puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
    puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
    puts $tcl_file "puts \$namd_file \"TMD on\""
    puts $tcl_file "puts \$namd_file \"TMDk $::comd::spring_k\""
    puts $tcl_file "puts \$namd_file \"TMDOutputFreq [expr $::comd::tmd_len*1]\""
    puts $tcl_file "puts \$namd_file \"TMDFile ..\/walker2_adjust.pdb\""
    puts $tcl_file "puts \$namd_file \"TMDFirstStep 0\""
    puts $tcl_file "puts \$namd_file \"TMDLastStep [expr $::comd::tmd_len*5]\""
    puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
    puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
    puts $tcl_file "puts \$namd_file \"outputEnergies [expr $::comd::tmd_len*1]\""
    puts $tcl_file "puts \$namd_file \"outputPressure [expr $::comd::tmd_len*1]\""
    puts $tcl_file "puts \$namd_file \"restartfreq [expr $::comd::tmd_len*1]\""
    puts $tcl_file "puts \$namd_file \"dcdfreq [expr $::comd::tmd_len*1]\""
    puts $tcl_file "puts \$namd_file \"xstfreq [expr $::comd::tmd_len*1]\""
    puts $tcl_file "puts \$namd_file \"run [expr $::comd::tmd_len*5]\""
    puts $tcl_file "close \$namd_file"
    puts $tcl_file "puts \$sh_file \"cd ${::comd::output_prefix}_walker2_pro\""
    puts $tcl_file "puts \$sh_file \"\\\$NAMD \+devices $::comd::gpus_selection2 \+ppn $processes_per_run pro.conf > pro\$\{cycle\}.log \&\""
    puts $tcl_file "puts \$sh_file \"cd ..\""
  }
  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "puts \"Now running TMD \$\{cycle\}\""
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"

  puts $tcl_file "if {\$status} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: \$output\""
  puts $tcl_file "exit"
  puts $tcl_file "}"

  puts $tcl_file "puts \"Finished TMD \$\{cycle\}\""

  # If any files are missing continue to the end of the cycle and end up going back to this cycle again
  puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker1_pro/walker1_process\$\{cycle\}.coor r} fid\]} {"
  puts $tcl_file "continue"
  puts $tcl_file "}"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker1_pro/walker1_process\$\{cycle\}.coor r} fid\]} {"
    puts $tcl_file "continue"
    puts $tcl_file "}"
  }

  ###### MINIMIZATION AT THE END OF THE LOOP ########
  puts $tcl_file "set sh_file \[open \"${::comd::output_prefix}_min_\$cycle.sh\" w\]"
  puts $tcl_file "set sh_filename \"${::comd::output_prefix}_min_\$cycle.sh\""
  puts $tcl_file "puts \$sh_file \"\\\#\\\!\\\/bin\\\/bash\""
  puts $tcl_file "puts \$sh_file \"NAMD=\\\"\$namd2path \\\"\""

  # Walker 1 minimization
  puts $tcl_file "set namd_file \[open \[file join \"${::comd::output_prefix}_walker1_min\" \"min.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ../walker1_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ../walker1_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $::comd::temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname walker1_minimized\${cycle}\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  
  foreach tempfile $::comd::para_file {
    puts $tcl_file "puts \$namd_file \"parameters $tempfile\""
  }
  
  puts $tcl_file "puts \$namd_file \"set restartname res\""
  puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${::comd::output_prefix}_walker1_pro\/walker1_process\$\{cycle\}.coor\""
  puts $tcl_file "puts \$namd_file \"binvelocities ..\/${::comd::output_prefix}_walker1_pro\/walker1_process\$\{cycle\}.vel\""
  puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${::comd::output_prefix}_walker1_pro\/walker1_process\$\{cycle\}.xst\""
  puts $tcl_file "puts \$namd_file \"wrapWater on\""
  puts $tcl_file "puts \$namd_file \"wrapAll on\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"PME yes\""
  puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies $::comd::min_length\""
  puts $tcl_file "puts \$namd_file \"outputPressure $::comd::min_length\""
  puts $tcl_file "puts \$namd_file \"restartfreq $::comd::min_length\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr ${::comd::min_length}*5]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr ${::comd::min_length}*5]\""
  puts $tcl_file "puts \$namd_file \"minimize [expr ${::comd::min_length}*5]\""
  puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${::comd::output_prefix}_walker1_min\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD \+devices $::comd::gpus_selection1 \+ppn $processes_per_run min.conf > min\$\{cycle\}.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\""

  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    # Walker 2 minimization
    puts $tcl_file "set namd_file \[open \[file join \"${::comd::output_prefix}_walker2_min\" \"min.conf\"\] w\]"
    puts $tcl_file "puts \$namd_file \"coordinates     ../walker2_ionized.pdb\""
    puts $tcl_file "puts \$namd_file \"structure       ../walker2_ionized.psf\""
    puts $tcl_file "puts \$namd_file \"set temperature $::comd::temperature\""
    puts $tcl_file "puts \$namd_file \"set outputname walker2_minimized\${cycle}\""
    puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
    puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
    
    foreach tempfile $::comd::para_file {
      puts $tcl_file "puts \$namd_file \"parameters $tempfile\""
    }

    puts $tcl_file "puts \$namd_file \"set restartname res\""
    puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${::comd::output_prefix}_walker2_pro\/walker2_process\$\{cycle\}.coor\""
    puts $tcl_file "puts \$namd_file \"binvelocities ..\/${::comd::output_prefix}_walker2_pro\/walker2_process\$\{cycle\}.vel\""
    puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${::comd::output_prefix}_walker2_pro\/walker2_process\$\{cycle\}.xst\""
    puts $tcl_file "puts \$namd_file \"wrapWater on\""
    puts $tcl_file "puts \$namd_file \"wrapAll on\""
    puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
    puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
    puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
    puts $tcl_file "puts \$namd_file \"timestep        1.0\""
    puts $tcl_file "puts \$namd_file \"switching       on\""
    puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
    puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
    puts $tcl_file "puts \$namd_file \"rigidBonds none\""
    puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
    puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
    puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
    puts $tcl_file "puts \$namd_file \"PME yes\""
    puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
    puts $tcl_file "puts \$namd_file \"langevin on\""
    puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
    puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
    puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
    puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
    puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
    puts $tcl_file "puts \$namd_file \"outputEnergies $::comd::min_length\""
    puts $tcl_file "puts \$namd_file \"outputPressure $::comd::min_length\""
    puts $tcl_file "puts \$namd_file \"restartfreq $::comd::min_length\""
    puts $tcl_file "puts \$namd_file \"dcdfreq [expr $::comd::min_length*5]\""
    puts $tcl_file "puts \$namd_file \"xstfreq [expr $::comd::min_length*5]\""
    puts $tcl_file "puts \$namd_file \"minimize [expr $::comd::min_length*5]\""
    puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
    puts $tcl_file "close \$namd_file"
    puts $tcl_file "puts \$sh_file \"cd ${::comd::output_prefix}_walker2_min\""
    puts $tcl_file "puts \$sh_file \"\\\$NAMD \+devices $::comd::gpus_selection2 \+ppn $processes_per_run min.conf > min\$\{cycle\}.log \&\""
    puts $tcl_file "puts \$sh_file \"cd ..\""
  }

  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "puts \"Now running minimization \$\{cycle\}\""
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"

  puts $tcl_file "if {\$status} {"
  puts $tcl_file "set err_file \[open \"$::comd::output_prefix.log\" a\]"
  puts $tcl_file "puts \$err_file \"ERROR: \$output\""
  puts $tcl_file "exit"
  puts $tcl_file "}"

  puts $tcl_file "puts \"Finished minimization \$\{cycle\}\""

  # Add the resulting PDBs to DCD files with the other ones from previous cycles
  puts $tcl_file "set status \[catch \{exec prody catdcd initr.dcd ${::comd::output_prefix}_walker1_min\/walker1_minimized\$\{cycle\}.dcd -o walker1_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv walker1_trajectory.dcd initr.dcd\} output\]" 
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "set status \[catch \{exec prody catdcd fintr.dcd ${::comd::output_prefix}_walker2_min\/walker2_minimized\$\{cycle\}.dcd -o walker2_trajectory.dcd\} output\]"
    puts $tcl_file "set status \[catch \{exec mv walker2_trajectory.dcd fintr.dcd\} output\]"
  }
  puts $tcl_file "puts \"Finished concatenating trajectories for cycle \$\{cycle\}\""

  # If files are missing continue to the end of the loop and the next loop will retry this cycle
  puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker1_min/walker1_minimized\$\{cycle\}.coor r} fid\]} {"
  puts $tcl_file "continue"
  puts $tcl_file "}"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    puts $tcl_file "if {\[catch {open ${::comd::output_prefix}_walker2_min/walker2_minimized\$\{cycle\}.coor r} fid\]} {"
    puts $tcl_file "continue"
    puts $tcl_file "}"
  }
  
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} {
    # calculate and output RMSD between the two endpoints if transitioning
    puts $tcl_file "mol delete all" 
    puts $tcl_file "mol load psf walker1_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker1_min/walker1_minimized\${cycle}.coor" 
    puts $tcl_file "set sel1 \[atomselect top \"name CA\"\]" 
    puts $tcl_file "set sel1a \[atomselect top all\]"
    puts $tcl_file "mol load psf walker2_ionized.psf"
    puts $tcl_file "mol addfile ${::comd::output_prefix}_walker2_min/walker2_minimized\${cycle}.coor"  
    puts $tcl_file "set sel2 \[atomselect top \"name CA\"\]" 
    puts $tcl_file "set sel2a \[atomselect top all\]"
    puts $tcl_file "set trans_mat \[measure fit \$sel2 \$sel1\]"
    puts $tcl_file "\$sel2a move \$trans_mat"
    puts $tcl_file "set rmsd \[measure rmsd \$sel2 \$sel1\]"
    puts $tcl_file "set all_rmsd(\$\{cycle\}) \$rmsd"
    puts $tcl_file "puts \$rmsd_file \"\$rmsd\""
    puts $tcl_file "if \{\(\$rmsd < 1.5)\|\|(\[expr \$all_rmsd\(\[expr \$\{cycle\}\-1\]\) - \$rmsd]\ < 0.15 \)\} \{ break \}"
  }

  # end loop
  puts $tcl_file "}"

  # tidy up
  puts $tcl_file "set status \[catch \{exec mv initr.dcd walker1_trajectory.dcd\} output\]"
  if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}]} { 
    puts $tcl_file "set status \[catch \{exec mv fintr.dcd walker2_trajectory.dcd\} output\]" 
  }
  #AJ
  puts $tcl_file "exit"
  close $tcl_file
  file delete pro_formatted.pdb
  file delete pro_formatted_autopsf.pdb
  file delete pro_formatted_autopsf.psf
  file delete pro_formatted_autopsf.log
  file delete pro_wb.psf
  file delete pro_wb.pdb
  file delete pro_wb.log
  file delete fino.pdb
  file delete init.pdb

  file copy -force $COMD_PATH/anmmc.py $::comd::outputdir/anmmc.py

  if {[info exists ::comd::from_commandline] == 0} {
    ::comd::Logview [file join "$::comd::outputdir" "$::comd::output_prefix.log"]
    tk_messageBox -type ok -title "Setup Complete" \
      -message "Setup of $::comd::output_prefix is complete. See $::comd::output_prefix.log file."
  }

  if {$::comd::run_now} {
    puts "Now running"
    cd $::comd::start_dir
    source $tcl_file_name
    puts "Finished"
  }

}

proc comd_tk {} {
  ::comd::comdgui
}

if { $argc < 3 } {
  puts "comd.tcl requires at least two arguments: filenames for the starting PDBs."
  puts "Please provide the same filename twice to calculate a random walk "
  puts "rather than a transition."
} else {

  if {[catch {
    set num_args 25

    # Take parameter values from input arguments as far as possible
    for {set index 0} {$index < $argc -1} {incr index} {
      if {$index eq  0} {set ::comd::outputdir [lindex $argv $index]}
      if {$index eq  1} {set ::comd::output_prefix [lindex $argv $index]}
      if {$index eq  2} {set ::comd::walker1_pdb [lindex $argv $index]}
      if {$index eq  3} {set ::comd::walker2_pdb [lindex $argv $index]}
      if {$index eq  4} {
        set ::comd::comd_cycle [lindex $argv $index]
	set ::comd::comd_cycle [expr ${::comd::comd_cycle}+1]
      }
      if {$index eq  5} {
        set ::comd::dev_mag [lindex $argv $index]
        set ::comd::dev_mag [expr $::comd::dev_mag]
      }
      if {$index eq  6} {
        set ::comd::accept_para [lindex $argv $index]
        set ::comd::accept_para [expr $::comd::accept_para]
      }
      if {$index eq  7} {
        set ::comd::step_cutoff [lindex $argv $index]
        set ::comd::step_cutoff [expr $::comd::step_cutoff]
      }
      if {$index eq  8} {
        set ::comd::min_length [lindex $argv $index]
        set ::comd::min_length [expr int($::comd::min_length * 100)]
      }
      if {$index eq  9} {
        set ::comd::tmd_len [lindex $argv $index]
        set ::comd::tmd_len [expr int($::comd::tmd_len * 100)]
      }
      if {$index eq 10} {set ::comd::anm_cutoff [lindex $argv $index]}
      if {$index eq 11} {set ::comd::max_steps [lindex $argv $index]}
      if {$index eq 12} {set ::comd::accept_para [lindex $argv $index]}
      if {$index eq 13} {set ::comd::walker1_chid [lindex $argv $index]}
      if {$index eq 14} {set ::comd::walker2_chid [lindex $argv $index]}
      if {$index eq 15} {set ::comd::solvent_padding_x [lindex $argv $index]}
      if {$index eq 16} {set ::comd::solvent_padding_y [lindex $argv $index]}
      if {$index eq 17} {set ::comd::solvent_padding_z [lindex $argv $index]}
      if {$index eq 18} {set ::comd::topo_file [lindex $argv $index]}
      if {$index eq 19} {set ::comd::temperature [lindex $argv $index]}
      if {$index eq 20} {set ::comd::para_file [list [lindex $argv $index]]}
      if {$index eq 21} {set ::comd::spring_k [lindex $argv $index]}
      if {$index eq 22} {
        set ::comd::gpus_selected [lindex $argv $index]
        set ::comd::gpus_present 1
      }
      if {$index eq 23} {set ::comd::num_cores [lindex $argv $index]}
      if {$index eq 24} {set ::comd::run_now [lindex $argv $index]}
    }

    # Fill in the remaining values with defaults
    for {set index $index} {$index < $num_args} {incr index} {
      if {$index eq  4} {set ::comd::comd_cycle 100}
      if {$index eq  5} {set ::comd::dev_mag 0}
      if {$index eq  6} {set ::comd::accept_para ""}
      if {$index eq  7} {set ::comd::step_cutoff 0}
      if {$index eq  8} {set ::comd::min_length 100}
      if {$index eq  9} {set ::comd::tmd_len 10}
      if {$index eq 10} {set ::comd::anm_cutoff ""}
      if {$index eq 11} {set ::comd::max_steps [lindex $argv $index]}
      if {$index eq 15} {set ::comd::solvent_padding_x 10}
      if {$index eq 16} {set ::comd::solvent_padding_y 10}
      if {$index eq 17} {set ::comd::solvent_padding_z 10}
      if {$index eq 18} {set ::comd::topo_file [list]}
      if {$index eq 19} {set ::comd::temperature 298}
      if {$index eq 20} {set ::comd::para_file [list]}
      if {$index eq 21} {set ::comd::spring_k 20000}
      if {$index eq 22} {
        if {[catch {
          set output [eval exec "nvidia-smi"]
          set records [split $output "\n"]

          set j [llength $records]

          set k 0
          set i 0
          set done_header 0
          set found_processes 0
          set ::comd::gpus_selected [list]
          foreach rec $records {

            set found_processes [string match *Processes* $rec]
            if {$found_processes == 1} {break}

            if {$i == 6 && $done_header == 0} {
              set done_header 1
              set i 0
            } elseif {$done_header && $i == 1} {
              set fields [split $rec " "]
              lappend ::comd::gpus_selected [lindex $fields 3]
            } elseif {$done_header && $i == 3} {
              set i 0
            }

            incr i
            incr k
          }

          set ::comd::gpus_selected [lreplace $::comd::gpus_selected [expr {[llength $::comd::gpus_selected]-1 }] [expr {[llength $::comd::gpus_selected]-1 }]]

          if {[expr {$::comd::walker1_pdb}] ne [expr {$::comd::walker2_pdb}] 
          && [expr [llength $::comd::gpus_selected] > 1]
          } then {
            set selection1 [list]
            set selection2 [list]
            for {set i 0} {$i < [expr [llength $::comd::gpus_selected]/2]} {incr i} {
              lappend selection1 [lindex $::comd::gpus_selected $i]
              lappend selection2 [lindex $::comd::gpus_selected [expr {${i} + [llength $::comd::gpus_selected]/2 }]]
            }
            set ::comd::gpus_selection1 [join $selection1 ","]
            set ::comd::gpus_selection2 [join $selection2 ","]
          } else {
            if {[expr [llength $::comd::gpus_selected] > 1]} {
              set ::comd::gpus_selected [join $::comd::gpus_selected ","]
            }
            set ::comd::gpus_selection1 $::comd::gpus_selected
            set ::comd::gpus_selection2 $::comd::gpus_selected
          }
        }]} {
          set ::comd::gpus_present 0
        } else {
          set ::comd::gpus_present 1
        }
      }
      if {$index eq 23} {
        set ::comd::num_cores [expr {[eval exec "cat /proc/cpuinfo | grep processor | tail -n 1 | awk \" \{ print \\\$3 \} \""] + 1}]}
      if {$index eq 24} {set ::comd::run_now 1}
      if {$index eq 25} {set ::comd::from_commandline 1}
    }

    set ::comd::start_dir [pwd]
    ::comd::Prepare_system
    exit

  } result ]} {
    puts "coMD simulation FINISHED successfully"
    set log_file [open [file join "$::comd::outputdir" "$::comd::output_prefix.log"] a]
    puts $log_file "coMD simulation FINISHED successfully"
    close $log_file
  } else {
    puts "$result"
    set log_file [open [file join "$::comd::outputdir" "$::comd::output_prefix.log"] a] 
    puts $log_file "ERROR: $result"
    close $log_file
  }
}

