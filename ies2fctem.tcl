#!/bin/sh
# the next line restarts using wish\
exec wish "$0" "$@" 

if {![info exists vTcl(sourcing)]} {

    # Provoke name search
    catch {package require bogus-package-name}
    set packageNames [package names]

    # BLT is needed
    package require BLT
    
    package require Tk
    switch $tcl_platform(platform) {
	windows {
            option add *Button.padY 0
	}
	default {
            option add *Scrollbar.width 10
            option add *Scrollbar.highlightThickness 0
            option add *Scrollbar.elementBorderWidth 2
            option add *Scrollbar.borderWidth 2
	}
    }
    
    # Tix is required
    package require Tix
    
}

#############################################################################
# Visual Tcl v1.60 Project
#




#############################################################################
## vTcl Code to Load Stock Images


if {![info exist vTcl(sourcing)]} {
#############################################################################
## Procedure:  vTcl:rename

proc ::vTcl:rename {name} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    regsub -all "\\." $name "_" ret
    regsub -all "\\-" $ret "_" ret
    regsub -all " " $ret "_" ret
    regsub -all "/" $ret "__" ret
    regsub -all "::" $ret "__" ret

    return [string tolower $ret]
}

#############################################################################
## Procedure:  vTcl:image:create_new_image

proc ::vTcl:image:create_new_image {filename {description {no description}} {type {}} {data {}}} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    # Does the image already exist?
    if {[info exists ::vTcl(images,files)]} {
        if {[lsearch -exact $::vTcl(images,files) $filename] > -1} { return }
    }

    if {![info exists ::vTcl(sourcing)] && [string length $data] > 0} {
        set object [image create  [vTcl:image:get_creation_type $filename]  -data $data]
    } else {
        # Wait a minute... Does the file actually exist?
        if {! [file exists $filename] } {
            # Try current directory
            set script [file dirname [info script]]
            set filename [file join $script [file tail $filename] ]
        }

        if {![file exists $filename]} {
            set description "file not found!"
            ## will add 'broken image' again when img is fixed, for now create empty
            set object [image create photo -width 1 -height 1]
        } else {
            set object [image create  [vTcl:image:get_creation_type $filename]  -file $filename]
        }
    }

    set reference [vTcl:rename $filename]
    set ::vTcl(images,$reference,image)       $object
    set ::vTcl(images,$reference,description) $description
    set ::vTcl(images,$reference,type)        $type
    set ::vTcl(images,filename,$object)       $filename

    lappend ::vTcl(images,files) $filename
    lappend ::vTcl(images,$type) $object

    # return image name in case caller might want it
    return $object
}

#############################################################################
## Procedure:  vTcl:image:get_image

proc ::vTcl:image:get_image {filename} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    set reference [vTcl:rename $filename]

    # Let's do some checking first
    if {![info exists ::vTcl(images,$reference,image)]} {
        # Well, the path may be wrong; in that case check
        # only the filename instead, without the path.

        set imageTail [file tail $filename]

        foreach oneFile $::vTcl(images,files) {
            if {[file tail $oneFile] == $imageTail} {
                set reference [vTcl:rename $oneFile]
                break
            }
        }
    }
    return $::vTcl(images,$reference,image)
}

#############################################################################
## Procedure:  vTcl:image:get_creation_type

proc ::vTcl:image:get_creation_type {filename} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    switch [string tolower [file extension $filename]] {
        .ppm -
        .jpg -
        .bmp -
        .gif    {return photo}
        .xbm    {return bitmap}
        default {return photo}
    }
}

foreach img {

        {{[file join / var lib vtcl-1.6.0 images edit remove.gif]} {} stock {}}

            } {
    eval set _file [lindex $img 0]
    vTcl:image:create_new_image\
        $_file [lindex $img 1] [lindex $img 2] [lindex $img 3]
}

}
#################################
# VTCL LIBRARY PROCEDURES
#

if {![info exists vTcl(sourcing)]} {
#############################################################################
## Library Procedure:  Window

proc ::Window {args} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    global vTcl
    foreach {cmd name newname} [lrange $args 0 2] {}
    set rest    [lrange $args 3 end]
    if {$name == "" || $cmd == ""} { return }
    if {$newname == ""} { set newname $name }
    if {$name == "."} { wm withdraw $name; return }
    set exists [winfo exists $newname]
    switch $cmd {
        show {
            if {$exists} {
                wm deiconify $newname
            } elseif {[info procs vTclWindow$name] != ""} {
                eval "vTclWindow$name $newname $rest"
            }
            if {[winfo exists $newname] && [wm state $newname] == "normal"} {
                vTcl:FireEvent $newname <<Show>>
            }
        }
        hide    {
            if {$exists} {
                wm withdraw $newname
                vTcl:FireEvent $newname <<Hide>>
                return}
        }
        iconify { if $exists {wm iconify $newname; return} }
        destroy { if $exists {destroy $newname; return} }
    }
}
#############################################################################
## Library Procedure:  vTcl:DefineAlias

proc ::vTcl:DefineAlias {target alias widgetProc top_or_alias cmdalias} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    global widget
    set widget($alias) $target
    set widget(rev,$target) $alias
    if {$cmdalias} {
        interp alias {} $alias {} $widgetProc $target
    }
    if {$top_or_alias != ""} {
        set widget($top_or_alias,$alias) $target
        if {$cmdalias} {
            interp alias {} $top_or_alias.$alias {} $widgetProc $target
        }
    }
}
#############################################################################
## Library Procedure:  vTcl:DoCmdOption

proc ::vTcl:DoCmdOption {target cmd} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    ## menus are considered toplevel windows
    set parent $target
    while {[winfo class $parent] == "Menu"} {
        set parent [winfo parent $parent]
    }

    regsub -all {\%widget} $cmd $target cmd
    regsub -all {\%top} $cmd [winfo toplevel $parent] cmd

    uplevel #0 [list eval $cmd]
}
#############################################################################
## Library Procedure:  vTcl:FireEvent

proc ::vTcl:FireEvent {target event {params {}}} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    ## The window may have disappeared
    if {![winfo exists $target]} return
    ## Process each binding tag, looking for the event
    foreach bindtag [bindtags $target] {
        set tag_events [bind $bindtag]
        set stop_processing 0
        foreach tag_event $tag_events {
            if {$tag_event == $event} {
                set bind_code [bind $bindtag $tag_event]
                foreach rep "\{%W $target\} $params" {
                    regsub -all [lindex $rep 0] $bind_code [lindex $rep 1] bind_code
                }
                set result [catch {uplevel #0 $bind_code} errortext]
                if {$result == 3} {
                    ## break exception, stop processing
                    set stop_processing 1
                } elseif {$result != 0} {
                    bgerror $errortext
                }
                break
            }
        }
        if {$stop_processing} {break}
    }
}
#############################################################################
## Library Procedure:  vTcl:Toplevel:WidgetProc

proc ::vTcl:Toplevel:WidgetProc {w args} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    if {[llength $args] == 0} {
        ## If no arguments, returns the path the alias points to
        return $w
    }
    set command [lindex $args 0]
    set args [lrange $args 1 end]
    switch -- [string tolower $command] {
        "setvar" {
            foreach {varname value} $args {}
            if {$value == ""} {
                return [set ::${w}::${varname}]
            } else {
                return [set ::${w}::${varname} $value]
            }
        }
        "hide" - "show" {
            Window [string tolower $command] $w
        }
        "showmodal" {
            ## modal dialog ends when window is destroyed
            Window show $w; raise $w
            grab $w; tkwait window $w; grab release $w
        }
        "startmodal" {
            ## ends when endmodal called
            Window show $w; raise $w
            set ::${w}::_modal 1
            grab $w; tkwait variable ::${w}::_modal; grab release $w
        }
        "endmodal" {
            ## ends modal dialog started with startmodal, argument is var name
            set ::${w}::_modal 0
            Window hide $w
        }
        default {
            uplevel $w $command $args
        }
    }
}
#############################################################################
## Library Procedure:  vTcl:WidgetProc

proc ::vTcl:WidgetProc {w args} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    if {[llength $args] == 0} {
        ## If no arguments, returns the path the alias points to
        return $w
    }

    set command [lindex $args 0]
    set args [lrange $args 1 end]
    uplevel $w $command $args
}
#############################################################################
## Library Procedure:  vTcl:toplevel

proc ::vTcl:toplevel {args} {
    ## This procedure may be used free of restrictions.
    ##    Exception added by Christian Gavin on 08/08/02.
    ## Other packages and widget toolkits have different licensing requirements.
    ##    Please read their license agreements for details.

    uplevel #0 eval toplevel $args
    set target [lindex $args 0]
    namespace eval ::$target {set _modal 0}
}
}


if {[info exists vTcl(sourcing)]} {

proc vTcl:project:info {} {
    set base .top66
    namespace eval ::widgets::$base {
        set set,origin 1
        set set,size 1
        set runvisible 1
    }
    namespace eval ::widgets::$base.tix69 {
        array set save {-background 1 -disabledforeground 1 -filebitmap 1 -highlightbackground 1 -label 1 -options 1}
    }
    namespace eval ::widgets::$base.sca70 {
        array set save {-bigincrement 1 -from 1 -label 1 -orient 1 -resolution 1 -tickinterval 1 -to 1 -variable 1}
    }
    namespace eval ::widgets::$base.lab77 {
        array set save {-foreground 1 -highlightcolor 1 -text 1}
    }
    namespace eval ::widgets::$base.lab78 {
        array set save {-foreground 1 -highlightcolor 1 -text 1}
    }
    namespace eval ::widgets::$base.lab80 {
        array set save {-foreground 1 -highlightcolor 1 -text 1}
    }
    set site_3_0 $base.lab80
    namespace eval ::widgets::$site_3_0.che81 {
        array set save {-offvalue 1 -onvalue 1 -text 1 -variable 1}
    }
    namespace eval ::widgets::$site_3_0.che82 {
        array set save {-text 1 -variable 1}
    }
    namespace eval ::widgets::$base.but83 {
        array set save {-command 1 -disabledforeground 1 -text 1}
    }
    namespace eval ::widgets::$base.lab89 {
        array set save {-foreground 1 -highlightcolor 1 -text 1}
    }
    namespace eval ::widgets::$base.tix70 {
        array set save {-highlightbackground 1}
        namespace eval subOptions {
            array set save {-anchor 1 -label 1}
        }
    }
    set site_5_page1 [$base.tix70 subwidget [lindex [$base.tix70 pages] 0]]
    namespace eval ::widgets::$site_5_page1 {
        array set save {-height 1 -highlightcolor 1 -visual 1 -width 1}
    }
    set site_5_0 $site_5_page1
    namespace eval ::widgets::$site_5_0.lab77 {
        array set save {-disabledforeground 1 -relief 1 -text 1 -textvariable 1}
    }
    namespace eval ::widgets::$site_5_0.cpd78 {
        array set save {-disabledforeground 1 -relief 1 -text 1 -textvariable 1}
    }
    namespace eval ::widgets::$site_5_0.cpd80 {
        array set save {-disabledforeground 1 -highlightbackground 1 -text 1}
    }
    namespace eval ::widgets::$site_5_0.cpd81 {
        array set save {-disabledforeground 1 -highlightbackground 1 -text 1}
    }
    set site_5_page2 [$base.tix70 subwidget [lindex [$base.tix70 pages] 1]]
    namespace eval ::widgets::$site_5_page2 {
        array set save {-height 1 -highlightcolor 1 -width 1}
    }
    set site_5_0 $site_5_page2
    namespace eval ::widgets::$site_5_0.lab82 {
        array set save {-disabledforeground 1 -image 1 -relief 1 -text 1}
    }
    set site_5_page3 [$base.tix70 subwidget [lindex [$base.tix70 pages] 2]]
    namespace eval ::widgets::$site_5_page3 {
        array set save {-height 1 -highlightcolor 1 -width 1}
    }
    set site_5_0 $site_5_page3
    namespace eval ::widgets::$site_5_0.lab73 {
        array set save {-disabledforeground 1 -relief 1 -text 1}
    }
    namespace eval ::widgets::$site_5_0.cpd74 {
        array set save {-disabledforeground 1 -relief 1 -text 1}
    }
    namespace eval ::widgets::$site_5_0.gra84 {
        array set save {-background 1 -barmode 1 -borderwidth 1 -data 1 -height 1 -plotpadx 1 -plotpady 1 -width 1}
    }
    namespace eval ::widgets::$site_5_0.lab66 {
        array set save {-disabledforeground 1 -text 1}
    }
    namespace eval ::widgets_bindings {
        set tagslist _TopLevel
    }
    namespace eval ::vTcl::modules::main {
        set procs {
            init
            main
        }
        set compounds {
        }
        set projectType single
    }
}
}

#################################
# USER DEFINED PROCEDURES
#
#############################################################################
## Procedure:  main

proc ::main {argc argv} {

}

#############################################################################
## Initialization Procedure:  init

proc ::init {argc argv} {

}
set comp 0
set prevangle 910
init $argc $argv

#################################
# VTCL GENERATED GUI PROCEDURES
#

proc vTclWindow. {base} {
    if {$base == ""} {
        set base .
    }

    ###################
    # CREATING WIDGETS
    ###################
    wm focusmodel $top passive
    wm geometry $top 1x1+0+0; update
    wm maxsize $top 1265 994
    wm minsize $top 1 1
    wm overrideredirect $top 0
    wm resizable $top 0 0
    wm withdraw $top
    wm title $top "vtcl.tcl #2"
    bindtags $top "$top Vtcl.tcl all"
    vTcl:FireEvent $top <<Create>>
    wm protocol $top WM_DELETE_WINDOW "vTcl:FireEvent $top <<DeleteWindow>>"

    ###################
    # SETTING GEOMETRY
    ###################

    vTcl:FireEvent $base <<Ready>>
}

proc vTclWindow.top66 {base} {
    if {$base == ""} {
        set base .top66
    }
    if {[winfo exists $base]} {
        wm deiconify $base; return
    }
    set top $base

    ###################
    # CREATING WIDGETS
    ###################
    vTcl:toplevel $top -class Toplevel \
        -highlightcolor black 
    wm focusmodel $top passive
    wm geometry $top 445x599+106+179; update
    wm maxsize $top 1265 994
    wm minsize $top 1 1
    wm overrideredirect $top 0
    wm resizable $top 0 0
    wm deiconify $top
    wm title $top "IES file tilt and integration system V0.6"
    vTcl:DefineAlias "$top" "Toplevel1" vTcl:Toplevel:WidgetProc "" 1
    bindtags $top "$top Toplevel all _TopLevel"
    vTcl:FireEvent $top <<Create>>
    wm protocol $top WM_DELETE_WINDOW "vTcl:FireEvent $top <<DeleteWindow>>"

    tixFileEntry $top.tix69 \
        -label {IES file:} -height 24 -options {label.anchor e} \
        -variable iesfile -width 416
    bind $top.tix69 <FocusIn> {
        focus %W.frame.entry
    }
    .top66.tix69 subwidget frame configure -highlightthickness 2
    scale $top.sca70 \
        -bigincrement 0.0 -from 0.0 -label {Tilt angle (deg)} \
        -orient horizontal -resolution 1.0 -tickinterval 0.0 -to 90.0 \
        -variable tiltangle 
    vTcl:DefineAlias "$top.sca70" "Scale1" vTcl:WidgetProc "Toplevel1" 1
    labelframe $top.lab77 \
        -foreground black -text Label -highlightcolor black 
    vTcl:DefineAlias "$top.lab77" "Labelframe1" vTcl:WidgetProc "Toplevel1" 1
    labelframe $top.lab78 \
        -foreground black -text Label -highlightcolor black 
    vTcl:DefineAlias "$top.lab78" "Labelframe2" vTcl:WidgetProc "Toplevel1" 1
    labelframe $top.lab80 \
        -foreground black -text {Luminaire style} -highlightcolor black 
    vTcl:DefineAlias "$top.lab80" "Labelframe4" vTcl:WidgetProc "Toplevel1" 1
    set site_3_0 $top.lab80
    checkbutton $site_3_0.che81 \
        -offvalue 1 -onvalue 0 -text Other -variable style -command {Scale1 configure -from -90.0 -to 90.0} \
        -selectcolor #ff0000
    vTcl:DefineAlias "$site_3_0.che81" "Checkbutton1" vTcl:WidgetProc "Toplevel1" 1
    checkbutton $site_3_0.che82 \
        -text Projector -variable style -command {Scale1 configure -from 0.0 -to 90.0} \
        -selectcolor #ff0000
    vTcl:DefineAlias "$site_3_0.che82" "Checkbutton2" vTcl:WidgetProc "Toplevel1" 1
    place $site_3_0.che81 \
        -in $site_3_0 -x 237 -y 18 -width 60 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che82 \
        -in $site_3_0 -x 72 -y 18 -width 80 -height 22 -anchor nw \
        -bordermode ignore 

    labelframe $top.lab89 \
        -foreground black -text Results -highlightcolor black 
    vTcl:DefineAlias "$top.lab89" "Labelframe3" vTcl:WidgetProc "Toplevel1" 1
    tixNoteBook $top.tix70 \
        -highlightbackground #d9d9d9 
    vTcl:DefineAlias "$top.tix70" "TixNoteBook1" vTcl:WidgetProc "Toplevel1" 1
    $top.tix70 add page1 \
        -anchor center -label {Integrated fluxes} 
    $top.tix70 add page2 \
        -anchor center -label {IES tilted function} 
    $top.tix70 add page3 \
        -anchor center -label {Horizontal average} 
    set site_5_page1 [$top.tix70 subwidget [lindex [$top.tix70 pages] 0]]
    label $site_5_page1.lab77 \
        -disabledforeground #a1a4a1 -relief sunken -text ? \
        -textvariable upflux -background #ffffff
    vTcl:DefineAlias "$site_5_page1.lab77" "Label11" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page1.lab77 \
        -in $site_5_page1 -x 180 -y 15 -height 15 -width 60 -anchor nw -bordermode ignore 
    label $site_5_page1.cpd78 \
        -disabledforeground #a1a4a1 -relief sunken -text ? \
        -textvariable downflux -background #ffffff
    vTcl:DefineAlias "$site_5_page1.cpd78" "Label12" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page1.cpd78 \
        -in $site_5_page1 -x 180 -y 50 -height 15 -width 60 -anchor nw -bordermode inside 
    label $site_5_page1.cpd80 \
        -disabledforeground #a1a4a1 -highlightbackground #d9d9d9 \
        -text {Upward flux fraction:} 
    vTcl:DefineAlias "$site_5_page1.cpd80" "Label10" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page1.cpd80 \
        -in $site_5_page1 -x 15 -y 15 -anchor nw -bordermode inside 
    label $site_5_page1.cpd81 \
        -disabledforeground #a1a4a1 -highlightbackground #d9d9d9 \
        -text {Downward flux fraction:} 
    vTcl:DefineAlias "$site_5_page1.cpd81" "Label13" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page1.cpd81 \
        -in $site_5_page1 -x 15 -y 50 -anchor nw -bordermode inside 
    set site_5_page2 [$top.tix70 subwidget [lindex [$top.tix70 pages] 1]]
    label $site_5_page2.lab82 \
        -disabledforeground #a1a4a1 -background #ffffff \
        -image [vTcl:image:get_image [file join / var lib vtcl-1.6.0 images edit remove.gif]] \
        -relief sunken -text label 
    vTcl:DefineAlias "$site_5_page2.lab82" "Label14" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab82 \
        -in $site_5_page2 -x 10 -y 10 -width 368 -height 320 -anchor nw \
        -bordermode ignore 
     label $site_5_page2.lab82.inside \
        -disabledforeground #a1a4a1 -background #fffca0 \
        -image [vTcl:image:get_image [file join / var lib vtcl-1.6.0 images edit remove.gif]] \
        -relief sunken -text label 
    vTcl:DefineAlias "$site_5_page2.lab82.inside" "Label14.inside" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab82.inside \
        -in $site_5_page2.lab82 -x 51 -y 24 -width 270 -height 270 -anchor nw \
        -bordermode ignore        
        
#
# IES tilted file image labels
#
    label $site_5_page2.lab83 \
        -disabledforeground #a1a4a1 \
        -relief groove -text East -background #ffffff
    vTcl:DefineAlias "$site_5_page2.lab83" "Label14a" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab83 \
        -in $site_5_page2 -x 305 -y 310 -width 30 -height 10 -anchor nw \
        -bordermode ignore         
    label $site_5_page2.lab84 \
        -disabledforeground #a1a4a1 \
        -relief groove -text West -background #ffffff
    vTcl:DefineAlias "$site_5_page2.lab84" "Label14b" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab84 \
        -in $site_5_page2 -x 60 -y 310 -width 30 -height 10 -anchor nw \
        -bordermode ignore         
     label $site_5_page2.lab85 \
        -disabledforeground #a1a4a1 \
        -relief groove -text South -background #ffffff
    vTcl:DefineAlias "$site_5_page2.lab85" "Label14c" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab85 \
        -in $site_5_page2 -x 183 -y 310 -width 30 -height 10 -anchor nw \
        -bordermode ignore         
      label $site_5_page2.lab86 \
        -disabledforeground #a1a4a1 \
        -relief groove -text Zenith -background #ffffff
    vTcl:DefineAlias "$site_5_page2.lab86" "Label14d" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab86 \
        -in $site_5_page2 -x 10 -y 30 -width 40 -height 10 -anchor nw \
        -bordermode ignore          
       label $site_5_page2.lab87 \
        -disabledforeground #a1a4a1 \
        -relief groove -text Nadir -background #ffffff
    vTcl:DefineAlias "$site_5_page2.lab87" "Label14e" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab87 \
        -in $site_5_page2 -x 10 -y 295 -width 40 -height 10 -anchor nw \
        -bordermode ignore  
        label $site_5_page2.lab88 \
        -disabledforeground #a1a4a1 \
        -relief groove -text Horizon -background #ffffff
    vTcl:DefineAlias "$site_5_page2.lab88" "Label14f" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page2.lab88 \
        -in $site_5_page2 -x 10 -y 162 -width 40 -height 10 -anchor nw \
        -bordermode ignore  
       
               
    set site_5_page3 [$top.tix70 subwidget [lindex [$top.tix70 pages] 2]]
    label $site_5_page3.lab73 \
        -disabledforeground #a1a4a1 -relief sunken -text label 
    vTcl:DefineAlias "$site_5_page3.lab73" "Label1" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page3.lab73 \
        -in $site_5_page3 -x -120 -y 260 -anchor nw -bordermode ignore 
    label $site_5_page3.cpd74 \
        -disabledforeground #a1a4a1 -relief sunken -text label 
    vTcl:DefineAlias "$site_5_page3.cpd74" "Label2" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page3.cpd74 \
        -in $site_5_page3 -x -110 -y 275 -anchor nw -bordermode inside 
    ::blt::graph $site_5_page3.gra84 \
        -barmode infront -borderwidth 0  \
        -height 320 -plotpadx {8 8} -plotpady {8 8} -width 385 -plotbackground #ffffff
    vTcl:DefineAlias "$site_5_page3.gra84" "Graph2" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page3.gra84 \
        -in $site_5_page3 -x 0 -y 0 -width 385 -height 295 -anchor nw \
        -bordermode ignore 
    label $site_5_page3.lab66 \
        -disabledforeground #a1a4a1 -text {Zenith angle (deg)} 
    vTcl:DefineAlias "$site_5_page3.lab66" "Label3" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page3.lab66 \
        -in $site_5_page3 -x 140 -y 295 -anchor nw -bordermode ignore 
#
# save the file as 
#       
    label $site_5_page3.lab67 \
        -disabledforeground #a1a4a1 -text {Save as file:} 
    vTcl:DefineAlias "$site_5_page3.lab67" "Label3a" vTcl:WidgetProc "Toplevel1" 1
    place $site_5_page3.lab67 \
        -in $site_5_page3 -x 10 -y 317 -anchor nw -bordermode ignore 
    entry $site_5_page3.e75 -background #ffffff -textvariable nomfichier           
    place $site_5_page3.e75 \
        -in $site_5_page3 -x 75 -y 315 -width 250 -anchor nw -bordermode ignore         
    button $site_5_page3.b75  -text Save \
        -command { if { $nomfichier != "" } {
                     if { $prevangle != "910" } {
                         exec mv fctem.txt $nomfichier
                     }   
                   } 
                 } 
                                 
    place $site_5_page3.b75 \
        -in $site_5_page3 -x 335 -y 314 -width 40 -height 20 -anchor nw -bordermode ignore

#
#  execute
#        
    button $top.but83 \
        -command { if { $iesfile != "" } {
                   set exten [file extension $iesfile]
                   if { $exten == ".ies" || $exten == ".IES" } {
                  
                      set folder [file dirname $iesfile]
                      cd $folder
                      puts "folder=$folder"
                      set ax ""
                      set ay ""
                      if  { $style == "1" } {
                          set tilt [expr -90+$tiltangle ] 
                          Label3 configure -text {Zenith angle (deg)}   
                          Scale1 configure -from 0.0 -to 90.0            
                      }
                      if  { $style == "0" } {
                           set tilt $tiltangle
                           Label3 configure -text {Zenith angle (deg)} 
                           Scale1 configure -from -90.0 -to 90.0              
                       }
                       exec ies2fctem.bash  $iesfile $tilt fctem.txt
                       exec convert iestilt.pgm iestilt.gif
                       exec rm -f iestilt.pgm
                       image create photo imagei -file [file join $folder/iestilt.gif] 
                       image create photo imagetilt
                       imagetilt copy imagei -zoom 3 3
                       Label14.inside configure -image imagetilt
                       exec rm -f iestilt.gif
   #  reading averaged light pattern               
                       set xyfile [open fctem.txt "r"]                    
                       set inEOF -1
                       while {[gets $xyfile xy] != $inEOF} {
                            foreach {y x} $xy {
                               lappend ax $x
                               lappend ay $y
                               puts "x= $x y=$y"
                            }  
                       }
                       close $xyfile  
   #  reading integrated fluxes
                       set fluxfile [open intflux.tmp "r"]                
                          gets $fluxfile upflux
                          gets $fluxfile downflux
                       close $fluxfile 
                       exec rm -f intflux.tmp               
                       if { $prevangle == "910" } {
                           Graph2 element create "<$tiltangle" -xdata $ax -ydata $ay 
                           Graph2 element configure "<$tiltangle" -symbol ""                      
                           set prevangle $tiltangle
                           set prev $tilt                        
                       }                                                        
                       if { $tilt != $prev } {
#                       set comp 0
                             Graph2 element delete "<$prevangle"   
                             Graph2 element create "<$tiltangle" -xdata $ax -ydata $ay  
                             Graph2 element configure "<$tiltangle" -symbol ""
                             set prevangle $tiltangle 
                             set prev $tilt                                                       
                        }
                       } } } \
        -disabledforeground #a1a4a1 -text {Execute / refresh}
    vTcl:DefineAlias "$top.but83" "Button1" vTcl:WidgetProc "Toplevel1" 1        
#    
#  label copyright
#    
    label $top.labcopy \
        -disabledforeground #a1a4a1 -relief sunken -text {GNU Public License, CopyRight Martin Aubé, MEMO Environnement 2005}
    place $top.labcopy \
        -in $top -x 17 -y 583 -anchor nw -bordermode ignore         
        
        
        
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.tix69 \
        -in $top -x 15 -y 10 -width 416 -height 24 -anchor nw \
        -bordermode ignore 
    place $top.sca70 \
        -in $top -x 10 -y 35 -width 421 -height 48 -anchor nw \
        -bordermode ignore 
    place $top.lab77 \
        -in $top -x 132 -y 344 -anchor nw -bordermode ignore 
    place $top.lab78 \
        -in $top -x 305 -y 345 -anchor nw -bordermode ignore 
    place $top.lab80 \
        -in $top -x 15 -y 100 -width 416 -height 46 -anchor nw \
        -bordermode ignore 
    place $top.but83 \
        -in $top -x 15 -y 155 -width 412 -height 26 -anchor nw \
        -bordermode ignore 
    place $top.lab89 \
        -in $top -x 15 -y 185 -width 416 -height 396 -anchor nw \
        -bordermode ignore 
    place $top.tix70 \
        -in $top -x 25 -y 200 -width 391 -height 371 -anchor nw \
        -bordermode ignore 

    vTcl:FireEvent $base <<Ready>>
}

#############################################################################
## Binding tag:  _TopLevel

bind "_TopLevel" <<Create>> {
    if {![info exists _topcount]} {set _topcount 0}; incr _topcount
}
bind "_TopLevel" <<DeleteWindow>> {
    if {[set ::%W::_modal]} {
                vTcl:Toplevel:WidgetProc %W endmodal
            } else {
                destroy %W; if {$_topcount == 0} {exit}
            }
}
bind "_TopLevel" <Destroy> {
    if {[winfo toplevel %W] == "%W"} {incr _topcount -1}
}

Window show .
Window show .top66

main $argc $argv
