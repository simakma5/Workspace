transcript on
if {[file exists rtl_work]} {
	vdel -lib rtl_work -all
}
vlib rtl_work
vmap work rtl_work

vcom -93 -work work {C:/Uni/MAM/hw02/project/digit_switch.vhd}
vcom -93 -work work {C:/Uni/MAM/hw02/project/rom.vhd}

