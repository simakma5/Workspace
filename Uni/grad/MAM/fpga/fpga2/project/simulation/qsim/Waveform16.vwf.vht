-- Copyright (C) 2017  Intel Corporation. All rights reserved.
-- Your use of Intel Corporation's design tools, logic functions 
-- and other software and tools, and its AMPP partner logic 
-- functions, and any output files from any of the foregoing 
-- (including device programming or simulation files), and any 
-- associated documentation or information are expressly subject 
-- to the terms and conditions of the Intel Program License 
-- Subscription Agreement, the Intel Quartus Prime License Agreement,
-- the Intel FPGA IP License Agreement, or other applicable license
-- agreement, including, without limitation, that your use is for
-- the sole purpose of programming logic devices manufactured by
-- Intel and sold by Intel or its authorized distributors.  Please
-- refer to the applicable agreement for further details.

-- *****************************************************************************
-- This file contains a Vhdl test bench with test vectors .The test vectors     
-- are exported from a vector file in the Quartus Waveform Editor and apply to  
-- the top level entity of the current Quartus project .The user can use this   
-- testbench to simulate his design using a third-party simulation tool .       
-- *****************************************************************************
-- Generated on "12/07/2022 17:15:45"
                                                             
-- Vhdl Test Bench(with test vectors) for design  :          hw02_simakma5
-- 
-- Simulation tool : 3rd Party
-- 

LIBRARY ieee;                                               
USE ieee.std_logic_1164.all;                                

ENTITY hw02_simakma5_vhd_vec_tst IS
END hw02_simakma5_vhd_vec_tst;
ARCHITECTURE hw02_simakma5_arch OF hw02_simakma5_vhd_vec_tst IS
-- constants                                                 
-- signals                                                   
SIGNAL A : STD_LOGIC;
SIGNAL B : STD_LOGIC;
SIGNAL C : STD_LOGIC;
SIGNAL CLK_IN : STD_LOGIC;
SIGNAL CLK_OUT : STD_LOGIC;
SIGNAL D : STD_LOGIC;
SIGNAL D0 : STD_LOGIC;
SIGNAL D1 : STD_LOGIC;
SIGNAL D2 : STD_LOGIC;
SIGNAL D3 : STD_LOGIC;
SIGNAL Data : STD_LOGIC;
SIGNAL h : STD_LOGIC;
SIGNAL KEY0 : STD_LOGIC;
SIGNAL KEY1 : STD_LOGIC;
SIGNAL LATCH : STD_LOGIC;
SIGNAL STROBE : STD_LOGIC;
SIGNAL SW0 : STD_LOGIC;
SIGNAL SW1 : STD_LOGIC;
SIGNAL test1us : STD_LOGIC;
SIGNAL test10us : STD_LOGIC;
SIGNAL test_bcd_CLK : STD_LOGIC;
SIGNAL test_CLK : STD_LOGIC;
SIGNAL test_Sync : STD_LOGIC;
COMPONENT hw02_simakma5
	PORT (
	A : OUT STD_LOGIC;
	B : OUT STD_LOGIC;
	C : OUT STD_LOGIC;
	CLK_IN : IN STD_LOGIC;
	CLK_OUT : OUT STD_LOGIC;
	D : OUT STD_LOGIC;
	D0 : OUT STD_LOGIC;
	D1 : OUT STD_LOGIC;
	D2 : OUT STD_LOGIC;
	D3 : OUT STD_LOGIC;
	Data : OUT STD_LOGIC;
	h : OUT STD_LOGIC;
	KEY0 : IN STD_LOGIC;
	KEY1 : IN STD_LOGIC;
	LATCH : OUT STD_LOGIC;
	STROBE : OUT STD_LOGIC;
	SW0 : IN STD_LOGIC;
	SW1 : IN STD_LOGIC;
	test1us : OUT STD_LOGIC;
	test10us : OUT STD_LOGIC;
	test_bcd_CLK : IN STD_LOGIC;
	test_CLK : IN STD_LOGIC;
	test_Sync : IN STD_LOGIC
	);
END COMPONENT;
BEGIN
	i1 : hw02_simakma5
	PORT MAP (
-- list connections between master ports and signals
	A => A,
	B => B,
	C => C,
	CLK_IN => CLK_IN,
	CLK_OUT => CLK_OUT,
	D => D,
	D0 => D0,
	D1 => D1,
	D2 => D2,
	D3 => D3,
	Data => Data,
	h => h,
	KEY0 => KEY0,
	KEY1 => KEY1,
	LATCH => LATCH,
	STROBE => STROBE,
	SW0 => SW0,
	SW1 => SW1,
	test1us => test1us,
	test10us => test10us,
	test_bcd_CLK => test_bcd_CLK,
	test_CLK => test_CLK,
	test_Sync => test_Sync
	);

-- CLK_IN
t_prcs_CLK_IN: PROCESS
BEGIN
LOOP
	CLK_IN <= '0';
	WAIT FOR 10000 ps;
	CLK_IN <= '1';
	WAIT FOR 10000 ps;
	IF (NOW >= 10000000 ps) THEN WAIT; END IF;
END LOOP;
END PROCESS t_prcs_CLK_IN;

-- KEY1
t_prcs_KEY1: PROCESS
BEGIN
	KEY1 <= '1';
WAIT;
END PROCESS t_prcs_KEY1;
END hw02_simakma5_arch;
