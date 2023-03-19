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

-- VENDOR "Altera"
-- PROGRAM "Quartus Prime"
-- VERSION "Version 17.1.0 Build 590 10/25/2017 SJ Lite Edition"

-- DATE "12/07/2022 19:19:52"

-- 
-- Device: Altera 10M50DAF484C7G Package FBGA484
-- 

-- 
-- This VHDL file should be used for ModelSim-Altera (VHDL) only
-- 

LIBRARY FIFTYFIVENM;
LIBRARY IEEE;
USE FIFTYFIVENM.FIFTYFIVENM_COMPONENTS.ALL;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY 	hard_block IS
    PORT (
	devoe : IN std_logic;
	devclrn : IN std_logic;
	devpor : IN std_logic
	);
END hard_block;

-- Design Ports Information
-- ~ALTERA_TMS~	=>  Location: PIN_H2,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_TCK~	=>  Location: PIN_G2,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_TDI~	=>  Location: PIN_L4,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_TDO~	=>  Location: PIN_M5,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- ~ALTERA_CONFIG_SEL~	=>  Location: PIN_H10,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- ~ALTERA_nCONFIG~	=>  Location: PIN_H9,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_nSTATUS~	=>  Location: PIN_G9,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_CONF_DONE~	=>  Location: PIN_F8,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default


ARCHITECTURE structure OF hard_block IS
SIGNAL gnd : std_logic := '0';
SIGNAL vcc : std_logic := '1';
SIGNAL unknown : std_logic := 'X';
SIGNAL ww_devoe : std_logic;
SIGNAL ww_devclrn : std_logic;
SIGNAL ww_devpor : std_logic;
SIGNAL \~ALTERA_TMS~~padout\ : std_logic;
SIGNAL \~ALTERA_TCK~~padout\ : std_logic;
SIGNAL \~ALTERA_TDI~~padout\ : std_logic;
SIGNAL \~ALTERA_CONFIG_SEL~~padout\ : std_logic;
SIGNAL \~ALTERA_nCONFIG~~padout\ : std_logic;
SIGNAL \~ALTERA_nSTATUS~~padout\ : std_logic;
SIGNAL \~ALTERA_CONF_DONE~~padout\ : std_logic;
SIGNAL \~ALTERA_TMS~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_TCK~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_TDI~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_CONFIG_SEL~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_nCONFIG~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_nSTATUS~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_CONF_DONE~~ibuf_o\ : std_logic;

BEGIN

ww_devoe <= devoe;
ww_devclrn <= devclrn;
ww_devpor <= devpor;
END structure;


LIBRARY ALTERA;
LIBRARY FIFTYFIVENM;
LIBRARY IEEE;
USE ALTERA.ALTERA_PRIMITIVES_COMPONENTS.ALL;
USE FIFTYFIVENM.FIFTYFIVENM_COMPONENTS.ALL;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY 	hw02_simakma5 IS
    PORT (
	STROBE : OUT std_logic;
	test_CLK : IN std_logic;
	test_Sync : IN std_logic;
	SW1 : IN std_logic;
	SW0 : IN std_logic;
	Data : OUT std_logic;
	D0 : OUT std_logic;
	D1 : OUT std_logic;
	D2 : OUT std_logic;
	D3 : OUT std_logic;
	test_bcd_CLK : IN std_logic;
	KEY1 : IN std_logic;
	CLK_OUT : OUT std_logic;
	LATCH : OUT std_logic;
	A : OUT std_logic;
	B : OUT std_logic;
	C : OUT std_logic;
	D : OUT std_logic;
	h : OUT std_logic;
	CLK_IN : IN std_logic;
	KEY0 : IN std_logic
	);
END hw02_simakma5;

-- Design Ports Information
-- STROBE	=>  Location: PIN_AB9,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- Data	=>  Location: PIN_AB17,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- D0	=>  Location: PIN_A8,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- D1	=>  Location: PIN_B11,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- D2	=>  Location: PIN_B12,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- D3	=>  Location: PIN_J11,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- KEY1	=>  Location: PIN_A7,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- CLK_OUT	=>  Location: PIN_AA12,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- LATCH	=>  Location: PIN_B10,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- A	=>  Location: PIN_A11,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- B	=>  Location: PIN_A9,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- C	=>  Location: PIN_H12,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- D	=>  Location: PIN_A10,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- h	=>  Location: PIN_C12,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- CLK_IN	=>  Location: PIN_P11,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- KEY0	=>  Location: PIN_B8,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- SW1	=>  Location: PIN_C11,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- SW0	=>  Location: PIN_C10,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- test_Sync	=>  Location: PIN_R11,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- test_bcd_CLK	=>  Location: PIN_C9,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- test_CLK	=>  Location: PIN_D12,	 I/O Standard: 2.5 V,	 Current Strength: Default


ARCHITECTURE structure OF hw02_simakma5 IS
SIGNAL gnd : std_logic := '0';
SIGNAL vcc : std_logic := '1';
SIGNAL unknown : std_logic := 'X';
SIGNAL devoe : std_logic := '1';
SIGNAL devclrn : std_logic := '1';
SIGNAL devpor : std_logic := '1';
SIGNAL ww_devoe : std_logic;
SIGNAL ww_devclrn : std_logic;
SIGNAL ww_devpor : std_logic;
SIGNAL ww_STROBE : std_logic;
SIGNAL ww_test_CLK : std_logic;
SIGNAL ww_test_Sync : std_logic;
SIGNAL ww_SW1 : std_logic;
SIGNAL ww_SW0 : std_logic;
SIGNAL ww_Data : std_logic;
SIGNAL ww_D0 : std_logic;
SIGNAL ww_D1 : std_logic;
SIGNAL ww_D2 : std_logic;
SIGNAL ww_D3 : std_logic;
SIGNAL ww_test_bcd_CLK : std_logic;
SIGNAL ww_KEY1 : std_logic;
SIGNAL ww_CLK_OUT : std_logic;
SIGNAL ww_LATCH : std_logic;
SIGNAL ww_A : std_logic;
SIGNAL ww_B : std_logic;
SIGNAL ww_C : std_logic;
SIGNAL ww_D : std_logic;
SIGNAL ww_h : std_logic;
SIGNAL ww_CLK_IN : std_logic;
SIGNAL ww_KEY0 : std_logic;
SIGNAL \~QUARTUS_CREATED_ADC1~_CHSEL_bus\ : std_logic_vector(4 DOWNTO 0);
SIGNAL \~QUARTUS_CREATED_ADC2~_CHSEL_bus\ : std_logic_vector(4 DOWNTO 0);
SIGNAL \inst5|inst5|7~clkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \inst5|inst3|7~clkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \inst8|inst~clkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \inst8|inst1|Mux1~clkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \inst5|inst4|7~clkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \test_Sync~inputclkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \inst5|inst6|7~clkctrl_INCLK_bus\ : std_logic_vector(3 DOWNTO 0);
SIGNAL \KEY1~input_o\ : std_logic;
SIGNAL \CLK_IN~input_o\ : std_logic;
SIGNAL \KEY0~input_o\ : std_logic;
SIGNAL \~QUARTUS_CREATED_GND~I_combout\ : std_logic;
SIGNAL \~QUARTUS_CREATED_UNVM~~busy\ : std_logic;
SIGNAL \~QUARTUS_CREATED_ADC1~~eoc\ : std_logic;
SIGNAL \~QUARTUS_CREATED_ADC2~~eoc\ : std_logic;
SIGNAL \STROBE~output_o\ : std_logic;
SIGNAL \Data~output_o\ : std_logic;
SIGNAL \D0~output_o\ : std_logic;
SIGNAL \D1~output_o\ : std_logic;
SIGNAL \D2~output_o\ : std_logic;
SIGNAL \D3~output_o\ : std_logic;
SIGNAL \CLK_OUT~output_o\ : std_logic;
SIGNAL \LATCH~output_o\ : std_logic;
SIGNAL \A~output_o\ : std_logic;
SIGNAL \B~output_o\ : std_logic;
SIGNAL \C~output_o\ : std_logic;
SIGNAL \D~output_o\ : std_logic;
SIGNAL \h~output_o\ : std_logic;
SIGNAL \test_CLK~input_o\ : std_logic;
SIGNAL \inst8|inst2|inst~0_combout\ : std_logic;
SIGNAL \test_Sync~input_o\ : std_logic;
SIGNAL \test_Sync~inputclkctrl_outclk\ : std_logic;
SIGNAL \inst8|inst2|inst~q\ : std_logic;
SIGNAL \inst8|inst2|inst1~0_combout\ : std_logic;
SIGNAL \inst8|inst2|inst1~feeder_combout\ : std_logic;
SIGNAL \inst8|inst2|inst1~q\ : std_logic;
SIGNAL \inst8|inst2|inst2~0_combout\ : std_logic;
SIGNAL \inst8|inst2|inst2~feeder_combout\ : std_logic;
SIGNAL \inst8|inst2|inst2~q\ : std_logic;
SIGNAL \inst8|inst2|inst9~combout\ : std_logic;
SIGNAL \inst8|inst2|inst4~0_combout\ : std_logic;
SIGNAL \inst8|inst2|inst4~feeder_combout\ : std_logic;
SIGNAL \inst8|inst2|inst4~q\ : std_logic;
SIGNAL \inst8|inst2|inst5~0_combout\ : std_logic;
SIGNAL \inst8|inst2|inst5~q\ : std_logic;
SIGNAL \inst8|inst~combout\ : std_logic;
SIGNAL \inst8|inst~clkctrl_outclk\ : std_logic;
SIGNAL \inst8|inst2|inst3~0_combout\ : std_logic;
SIGNAL \inst8|inst2|inst3~feeder_combout\ : std_logic;
SIGNAL \inst8|inst2|inst3~q\ : std_logic;
SIGNAL \inst8|inst1|Mux0~0_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux0~1_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux0~2_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux0~3_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~0_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux0~4_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~1_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~2_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~3_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~4_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~5_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~6_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~combout\ : std_logic;
SIGNAL \inst8|inst1|Mux1~clkctrl_outclk\ : std_logic;
SIGNAL \SW1~input_o\ : std_logic;
SIGNAL \inst8|inst4|state.digit3~q\ : std_logic;
SIGNAL \inst8|inst4|state.digit0~0_combout\ : std_logic;
SIGNAL \inst8|inst4|state.digit0~q\ : std_logic;
SIGNAL \inst8|inst4|state.digit1~0_combout\ : std_logic;
SIGNAL \inst8|inst4|state.digit1~q\ : std_logic;
SIGNAL \inst8|inst4|state.digit2~q\ : std_logic;
SIGNAL \SW0~input_o\ : std_logic;
SIGNAL \inst8|inst4|h~0_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux2~0_combout\ : std_logic;
SIGNAL \inst8|inst1|Mux2~1_combout\ : std_logic;
SIGNAL \test_bcd_CLK~input_o\ : std_logic;
SIGNAL \inst5|inst3|7~0_combout\ : std_logic;
SIGNAL \inst5|inst3|7~feeder_combout\ : std_logic;
SIGNAL \inst5|inst3|7~q\ : std_logic;
SIGNAL \inst5|inst3|7~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst3|31~combout\ : std_logic;
SIGNAL \inst5|inst3|11~q\ : std_logic;
SIGNAL \inst5|inst3|14~0_combout\ : std_logic;
SIGNAL \inst5|inst3|14~q\ : std_logic;
SIGNAL \inst5|inst3|30~combout\ : std_logic;
SIGNAL \inst5|inst3|19~q\ : std_logic;
SIGNAL \inst5|inst4|7~0_combout\ : std_logic;
SIGNAL \inst5|inst4|7~feeder_combout\ : std_logic;
SIGNAL \inst5|inst4|7~q\ : std_logic;
SIGNAL \inst5|inst4|7~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst4|31~combout\ : std_logic;
SIGNAL \inst5|inst4|11~feeder_combout\ : std_logic;
SIGNAL \inst5|inst4|11~q\ : std_logic;
SIGNAL \inst5|inst4|14~0_combout\ : std_logic;
SIGNAL \inst5|inst4|14~q\ : std_logic;
SIGNAL \inst5|inst4|30~combout\ : std_logic;
SIGNAL \inst5|inst4|19~q\ : std_logic;
SIGNAL \inst5|inst5|7~0_combout\ : std_logic;
SIGNAL \inst5|inst5|7~feeder_combout\ : std_logic;
SIGNAL \inst5|inst5|7~q\ : std_logic;
SIGNAL \inst5|inst5|7~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst5|31~combout\ : std_logic;
SIGNAL \inst5|inst5|11~feeder_combout\ : std_logic;
SIGNAL \inst5|inst5|11~q\ : std_logic;
SIGNAL \inst5|inst5|14~0_combout\ : std_logic;
SIGNAL \inst5|inst5|14~q\ : std_logic;
SIGNAL \inst5|inst5|30~combout\ : std_logic;
SIGNAL \inst5|inst5|19~q\ : std_logic;
SIGNAL \inst5|inst6|7~0_combout\ : std_logic;
SIGNAL \inst5|inst6|7~feeder_combout\ : std_logic;
SIGNAL \inst5|inst6|7~q\ : std_logic;
SIGNAL \inst5|inst6|7~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst6|31~combout\ : std_logic;
SIGNAL \inst5|inst6|11~q\ : std_logic;
SIGNAL \inst5|inst6|14~0_combout\ : std_logic;
SIGNAL \inst5|inst6|14~q\ : std_logic;
SIGNAL \inst5|inst6|30~combout\ : std_logic;
SIGNAL \inst5|inst6|19~q\ : std_logic;
SIGNAL \inst5|inst13~1_combout\ : std_logic;
SIGNAL \inst5|inst13~0_combout\ : std_logic;
SIGNAL \inst5|inst13~combout\ : std_logic;
SIGNAL \inst5|inst1~0_combout\ : std_logic;
SIGNAL \inst5|inst1~1_combout\ : std_logic;
SIGNAL \inst5|inst1~combout\ : std_logic;
SIGNAL \inst5|inst~0_combout\ : std_logic;
SIGNAL \inst5|inst~1_combout\ : std_logic;
SIGNAL \inst5|inst~combout\ : std_logic;
SIGNAL \inst5|inst12~0_combout\ : std_logic;
SIGNAL \inst5|inst12~1_combout\ : std_logic;
SIGNAL \inst5|inst12~combout\ : std_logic;
SIGNAL \inst7|106~0_combout\ : std_logic;
SIGNAL \inst6|84~0_combout\ : std_logic;
SIGNAL \inst6|84~1_combout\ : std_logic;
SIGNAL \inst7|103~0_combout\ : std_logic;
SIGNAL \inst7|79~0_combout\ : std_logic;
SIGNAL \inst7|78~0_combout\ : std_logic;
SIGNAL \inst7|17~combout\ : std_logic;
SIGNAL \inst7|10~q\ : std_logic;
SIGNAL \inst7|22~0_combout\ : std_logic;
SIGNAL \inst7|23~q\ : std_logic;
SIGNAL \inst7|25~0_combout\ : std_logic;
SIGNAL \inst7|24~q\ : std_logic;
SIGNAL \inst7|45~0_combout\ : std_logic;
SIGNAL \inst7|47~q\ : std_logic;
SIGNAL \inst7|44~0_combout\ : std_logic;
SIGNAL \inst7|46~q\ : std_logic;
SIGNAL \inst7|50~0_combout\ : std_logic;
SIGNAL \inst7|48~q\ : std_logic;
SIGNAL \inst7|51~0_combout\ : std_logic;
SIGNAL \inst7|49~q\ : std_logic;
SIGNAL \inst7|73~0_combout\ : std_logic;
SIGNAL \inst7|75~q\ : std_logic;
SIGNAL \inst7|72~0_combout\ : std_logic;
SIGNAL \inst7|72~1_combout\ : std_logic;
SIGNAL \inst7|74~q\ : std_logic;
SIGNAL \inst7|78~1_combout\ : std_logic;
SIGNAL \inst7|76~q\ : std_logic;
SIGNAL \inst7|79~1_combout\ : std_logic;
SIGNAL \inst7|77~q\ : std_logic;
SIGNAL \inst7|103~q\ : std_logic;
SIGNAL \inst7|100~0_combout\ : std_logic;
SIGNAL \inst7|100~1_combout\ : std_logic;
SIGNAL \inst7|102~q\ : std_logic;
SIGNAL \inst7|106~1_combout\ : std_logic;
SIGNAL \inst7|104~q\ : std_logic;
SIGNAL \inst7|107~0_combout\ : std_logic;
SIGNAL \inst7|107~1_combout\ : std_logic;
SIGNAL \inst7|105~q\ : std_logic;
SIGNAL \inst7|123~0_combout\ : std_logic;
SIGNAL \inst7|124~q\ : std_logic;
SIGNAL \inst101~combout\ : std_logic;
SIGNAL \inst5|inst6|ALT_INV_7~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst4|ALT_INV_7~clkctrl_outclk\ : std_logic;
SIGNAL \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst3|ALT_INV_7~clkctrl_outclk\ : std_logic;
SIGNAL \inst5|inst5|ALT_INV_7~clkctrl_outclk\ : std_logic;
SIGNAL \ALT_INV_test_bcd_CLK~input_o\ : std_logic;
SIGNAL \inst5|inst4|ALT_INV_19~q\ : std_logic;
SIGNAL \inst5|inst3|ALT_INV_19~q\ : std_logic;
SIGNAL \inst5|inst5|ALT_INV_19~q\ : std_logic;
SIGNAL \inst8|inst1|ALT_INV_Mux2~1_combout\ : std_logic;
SIGNAL \inst8|inst4|ALT_INV_state.digit0~q\ : std_logic;

COMPONENT hard_block
    PORT (
	devoe : IN std_logic;
	devclrn : IN std_logic;
	devpor : IN std_logic);
END COMPONENT;

BEGIN

STROBE <= ww_STROBE;
ww_test_CLK <= test_CLK;
ww_test_Sync <= test_Sync;
ww_SW1 <= SW1;
ww_SW0 <= SW0;
Data <= ww_Data;
D0 <= ww_D0;
D1 <= ww_D1;
D2 <= ww_D2;
D3 <= ww_D3;
ww_test_bcd_CLK <= test_bcd_CLK;
ww_KEY1 <= KEY1;
CLK_OUT <= ww_CLK_OUT;
LATCH <= ww_LATCH;
A <= ww_A;
B <= ww_B;
C <= ww_C;
D <= ww_D;
h <= ww_h;
ww_CLK_IN <= CLK_IN;
ww_KEY0 <= KEY0;
ww_devoe <= devoe;
ww_devclrn <= devclrn;
ww_devpor <= devpor;

\~QUARTUS_CREATED_ADC1~_CHSEL_bus\ <= (\~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\);

\~QUARTUS_CREATED_ADC2~_CHSEL_bus\ <= (\~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\);

\inst5|inst5|7~clkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \inst5|inst5|7~q\);

\inst5|inst3|7~clkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \inst5|inst3|7~q\);

\inst8|inst~clkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \inst8|inst~combout\);

\inst8|inst1|Mux1~clkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \inst8|inst1|Mux1~combout\);

\inst5|inst4|7~clkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \inst5|inst4|7~q\);

\test_Sync~inputclkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \test_Sync~input_o\);

\inst5|inst6|7~clkctrl_INCLK_bus\ <= (vcc & vcc & vcc & \inst5|inst6|7~q\);
\inst5|inst6|ALT_INV_7~clkctrl_outclk\ <= NOT \inst5|inst6|7~clkctrl_outclk\;
\inst5|inst4|ALT_INV_7~clkctrl_outclk\ <= NOT \inst5|inst4|7~clkctrl_outclk\;
\inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\ <= NOT \inst8|inst1|Mux1~clkctrl_outclk\;
\inst5|inst3|ALT_INV_7~clkctrl_outclk\ <= NOT \inst5|inst3|7~clkctrl_outclk\;
\inst5|inst5|ALT_INV_7~clkctrl_outclk\ <= NOT \inst5|inst5|7~clkctrl_outclk\;
\ALT_INV_test_bcd_CLK~input_o\ <= NOT \test_bcd_CLK~input_o\;
\inst5|inst4|ALT_INV_19~q\ <= NOT \inst5|inst4|19~q\;
\inst5|inst3|ALT_INV_19~q\ <= NOT \inst5|inst3|19~q\;
\inst5|inst5|ALT_INV_19~q\ <= NOT \inst5|inst5|19~q\;
\inst8|inst1|ALT_INV_Mux2~1_combout\ <= NOT \inst8|inst1|Mux2~1_combout\;
\inst8|inst4|ALT_INV_state.digit0~q\ <= NOT \inst8|inst4|state.digit0~q\;
auto_generated_inst : hard_block
PORT MAP (
	devoe => ww_devoe,
	devclrn => ww_devclrn,
	devpor => ww_devpor);

-- Location: LCCOMB_X44_Y52_N18
\~QUARTUS_CREATED_GND~I\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \~QUARTUS_CREATED_GND~I_combout\ = GND

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	combout => \~QUARTUS_CREATED_GND~I_combout\);

-- Location: IOOBUF_X34_Y0_N16
\STROBE~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst1|Mux0~3_combout\,
	devoe => ww_devoe,
	o => \STROBE~output_o\);

-- Location: IOOBUF_X69_Y0_N23
\Data~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst7|124~q\,
	devoe => ww_devoe,
	o => \Data~output_o\);

-- Location: IOOBUF_X46_Y54_N2
\D0~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst4|ALT_INV_state.digit0~q\,
	devoe => ww_devoe,
	o => \D0~output_o\);

-- Location: IOOBUF_X49_Y54_N9
\D1~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst4|state.digit1~q\,
	devoe => ww_devoe,
	o => \D1~output_o\);

-- Location: IOOBUF_X49_Y54_N2
\D2~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst4|state.digit2~q\,
	devoe => ww_devoe,
	o => \D2~output_o\);

-- Location: IOOBUF_X49_Y54_N23
\D3~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst4|state.digit3~q\,
	devoe => ww_devoe,
	o => \D3~output_o\);

-- Location: IOOBUF_X40_Y0_N2
\CLK_OUT~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst101~combout\,
	devoe => ww_devoe,
	o => \CLK_OUT~output_o\);

-- Location: IOOBUF_X46_Y54_N9
\LATCH~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst1|Mux2~1_combout\,
	devoe => ww_devoe,
	o => \LATCH~output_o\);

-- Location: IOOBUF_X51_Y54_N9
\A~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst5|inst1~combout\,
	devoe => ww_devoe,
	o => \A~output_o\);

-- Location: IOOBUF_X46_Y54_N23
\B~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst5|inst12~combout\,
	devoe => ww_devoe,
	o => \B~output_o\);

-- Location: IOOBUF_X49_Y54_N16
\C~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst5|inst~combout\,
	devoe => ww_devoe,
	o => \C~output_o\);

-- Location: IOOBUF_X51_Y54_N16
\D~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst5|inst13~combout\,
	devoe => ww_devoe,
	o => \D~output_o\);

-- Location: IOOBUF_X54_Y54_N30
\h~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \inst8|inst4|h~0_combout\,
	devoe => ww_devoe,
	o => \h~output_o\);

-- Location: IOIBUF_X51_Y54_N1
\test_CLK~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_test_CLK,
	o => \test_CLK~input_o\);

-- Location: LCCOMB_X47_Y48_N14
\inst8|inst2|inst~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst~0_combout\ = !\inst8|inst2|inst~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111100001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst8|inst2|inst~q\,
	combout => \inst8|inst2|inst~0_combout\);

-- Location: IOIBUF_X31_Y0_N1
\test_Sync~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_test_Sync,
	o => \test_Sync~input_o\);

-- Location: CLKCTRL_G18
\test_Sync~inputclkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \test_Sync~inputclkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \test_Sync~inputclkctrl_outclk\);

-- Location: FF_X47_Y48_N15
\inst8|inst2|inst\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst~clkctrl_outclk\,
	d => \inst8|inst2|inst~0_combout\,
	clrn => \test_Sync~inputclkctrl_outclk\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst2|inst~q\);

-- Location: LCCOMB_X46_Y48_N30
\inst8|inst2|inst1~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst1~0_combout\ = \inst8|inst2|inst~q\ $ (\inst8|inst2|inst1~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111111110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst8|inst2|inst~q\,
	datad => \inst8|inst2|inst1~q\,
	combout => \inst8|inst2|inst1~0_combout\);

-- Location: LCCOMB_X46_Y48_N10
\inst8|inst2|inst1~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst1~feeder_combout\ = \inst8|inst2|inst1~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst1~0_combout\,
	combout => \inst8|inst2|inst1~feeder_combout\);

-- Location: FF_X46_Y48_N11
\inst8|inst2|inst1\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst~clkctrl_outclk\,
	d => \inst8|inst2|inst1~feeder_combout\,
	clrn => \test_Sync~inputclkctrl_outclk\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst2|inst1~q\);

-- Location: LCCOMB_X46_Y48_N4
\inst8|inst2|inst2~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst2~0_combout\ = \inst8|inst2|inst2~q\ $ (((\inst8|inst2|inst~q\ & \inst8|inst2|inst1~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0011110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst2~q\,
	datac => \inst8|inst2|inst~q\,
	datad => \inst8|inst2|inst1~q\,
	combout => \inst8|inst2|inst2~0_combout\);

-- Location: LCCOMB_X46_Y48_N16
\inst8|inst2|inst2~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst2~feeder_combout\ = \inst8|inst2|inst2~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1010101010101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst2~0_combout\,
	combout => \inst8|inst2|inst2~feeder_combout\);

-- Location: FF_X46_Y48_N17
\inst8|inst2|inst2\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst~clkctrl_outclk\,
	d => \inst8|inst2|inst2~feeder_combout\,
	clrn => \test_Sync~inputclkctrl_outclk\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst2|inst2~q\);

-- Location: LCCOMB_X46_Y48_N18
\inst8|inst2|inst9\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst9~combout\ = (\inst8|inst2|inst1~q\ & (\inst8|inst2|inst2~q\ & (\inst8|inst2|inst~q\ & \inst8|inst2|inst3~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1000000000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datab => \inst8|inst2|inst2~q\,
	datac => \inst8|inst2|inst~q\,
	datad => \inst8|inst2|inst3~q\,
	combout => \inst8|inst2|inst9~combout\);

-- Location: LCCOMB_X46_Y48_N14
\inst8|inst2|inst4~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst4~0_combout\ = \inst8|inst2|inst9~combout\ $ (\inst8|inst2|inst4~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111111110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst8|inst2|inst9~combout\,
	datad => \inst8|inst2|inst4~q\,
	combout => \inst8|inst2|inst4~0_combout\);

-- Location: LCCOMB_X46_Y48_N8
\inst8|inst2|inst4~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst4~feeder_combout\ = \inst8|inst2|inst4~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000011110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst8|inst2|inst4~0_combout\,
	combout => \inst8|inst2|inst4~feeder_combout\);

-- Location: FF_X46_Y48_N9
\inst8|inst2|inst4\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst~clkctrl_outclk\,
	d => \inst8|inst2|inst4~feeder_combout\,
	clrn => \test_Sync~inputclkctrl_outclk\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst2|inst4~q\);

-- Location: LCCOMB_X46_Y48_N26
\inst8|inst2|inst5~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst5~0_combout\ = \inst8|inst2|inst5~q\ $ (((\inst8|inst2|inst4~q\ & \inst8|inst2|inst9~combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0110011011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst4~q\,
	datab => \inst8|inst2|inst5~q\,
	datad => \inst8|inst2|inst9~combout\,
	combout => \inst8|inst2|inst5~0_combout\);

-- Location: FF_X46_Y48_N29
\inst8|inst2|inst5\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst~combout\,
	asdata => \inst8|inst2|inst5~0_combout\,
	clrn => \test_Sync~inputclkctrl_outclk\,
	sload => VCC,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst2|inst5~q\);

-- Location: LCCOMB_X47_Y48_N24
\inst8|inst\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst~combout\ = LCELL((\test_CLK~input_o\ & ((\inst8|inst2|inst4~q\) # ((!\inst8|inst1|Mux0~2_combout\) # (!\inst8|inst2|inst5~q\)))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1000101010101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \test_CLK~input_o\,
	datab => \inst8|inst2|inst4~q\,
	datac => \inst8|inst2|inst5~q\,
	datad => \inst8|inst1|Mux0~2_combout\,
	combout => \inst8|inst~combout\);

-- Location: CLKCTRL_G12
\inst8|inst~clkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \inst8|inst~clkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \inst8|inst~clkctrl_outclk\);

-- Location: LCCOMB_X46_Y48_N12
\inst8|inst2|inst3~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst3~0_combout\ = \inst8|inst2|inst3~q\ $ (((\inst8|inst2|inst1~q\ & (\inst8|inst2|inst2~q\ & \inst8|inst2|inst~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0111111110000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datab => \inst8|inst2|inst2~q\,
	datac => \inst8|inst2|inst~q\,
	datad => \inst8|inst2|inst3~q\,
	combout => \inst8|inst2|inst3~0_combout\);

-- Location: LCCOMB_X46_Y48_N6
\inst8|inst2|inst3~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst2|inst3~feeder_combout\ = \inst8|inst2|inst3~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst3~0_combout\,
	combout => \inst8|inst2|inst3~feeder_combout\);

-- Location: FF_X46_Y48_N7
\inst8|inst2|inst3\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst~clkctrl_outclk\,
	d => \inst8|inst2|inst3~feeder_combout\,
	clrn => \test_Sync~inputclkctrl_outclk\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst2|inst3~q\);

-- Location: LCCOMB_X47_Y48_N26
\inst8|inst1|Mux0~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux0~0_combout\ = (\inst8|inst2|inst1~q\ & \inst8|inst2|inst2~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1010000010100000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datac => \inst8|inst2|inst2~q\,
	combout => \inst8|inst1|Mux0~0_combout\);

-- Location: LCCOMB_X47_Y48_N28
\inst8|inst1|Mux0~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux0~1_combout\ = (!\inst8|inst2|inst2~q\ & ((!\inst8|inst2|inst~q\) # (!\inst8|inst2|inst1~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000010101010101",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst2~q\,
	datac => \inst8|inst2|inst1~q\,
	datad => \inst8|inst2|inst~q\,
	combout => \inst8|inst1|Mux0~1_combout\);

-- Location: LCCOMB_X47_Y48_N6
\inst8|inst1|Mux0~2\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux0~2_combout\ = (\inst8|inst2|inst3~q\ & ((\inst8|inst1|Mux0~1_combout\))) # (!\inst8|inst2|inst3~q\ & (\inst8|inst1|Mux0~0_combout\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111110000110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst3~q\,
	datac => \inst8|inst1|Mux0~0_combout\,
	datad => \inst8|inst1|Mux0~1_combout\,
	combout => \inst8|inst1|Mux0~2_combout\);

-- Location: LCCOMB_X46_Y48_N28
\inst8|inst1|Mux0~3\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux0~3_combout\ = (\inst8|inst1|Mux0~2_combout\ & (\inst8|inst2|inst5~q\ & !\inst8|inst2|inst4~q\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000011000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst1|Mux0~2_combout\,
	datac => \inst8|inst2|inst5~q\,
	datad => \inst8|inst2|inst4~q\,
	combout => \inst8|inst1|Mux0~3_combout\);

-- Location: LCCOMB_X47_Y48_N8
\inst8|inst1|Mux1~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~0_combout\ = (\inst8|inst2|inst~q\) # (!\inst8|inst2|inst1~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111010111110101",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datac => \inst8|inst2|inst~q\,
	combout => \inst8|inst1|Mux1~0_combout\);

-- Location: LCCOMB_X47_Y48_N0
\inst8|inst1|Mux0~4\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux0~4_combout\ = (!\inst8|inst2|inst~q\) # (!\inst8|inst2|inst1~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0101111101011111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datac => \inst8|inst2|inst~q\,
	combout => \inst8|inst1|Mux0~4_combout\);

-- Location: LCCOMB_X47_Y48_N20
\inst8|inst1|Mux1~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~1_combout\ = (\inst8|inst2|inst2~q\ & (!\inst8|inst1|Mux1~0_combout\)) # (!\inst8|inst2|inst2~q\ & ((!\inst8|inst1|Mux0~4_combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000101001011111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst2~q\,
	datac => \inst8|inst1|Mux1~0_combout\,
	datad => \inst8|inst1|Mux0~4_combout\,
	combout => \inst8|inst1|Mux1~1_combout\);

-- Location: LCCOMB_X47_Y48_N30
\inst8|inst1|Mux1~2\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~2_combout\ = (\inst8|inst2|inst3~q\ & (\inst8|inst2|inst~q\)) # (!\inst8|inst2|inst3~q\ & ((!\inst8|inst1|Mux1~1_combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100000011110011",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst3~q\,
	datac => \inst8|inst2|inst~q\,
	datad => \inst8|inst1|Mux1~1_combout\,
	combout => \inst8|inst1|Mux1~2_combout\);

-- Location: LCCOMB_X47_Y48_N4
\inst8|inst1|Mux1~3\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~3_combout\ = (\inst8|inst2|inst4~q\ & ((\inst8|inst2|inst~q\))) # (!\inst8|inst2|inst4~q\ & (\inst8|inst1|Mux1~2_combout\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111110000110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst4~q\,
	datac => \inst8|inst1|Mux1~2_combout\,
	datad => \inst8|inst2|inst~q\,
	combout => \inst8|inst1|Mux1~3_combout\);

-- Location: LCCOMB_X47_Y48_N22
\inst8|inst1|Mux1~4\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~4_combout\ = (\inst8|inst2|inst1~q\) # (\inst8|inst2|inst~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111101011111010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datac => \inst8|inst2|inst~q\,
	combout => \inst8|inst1|Mux1~4_combout\);

-- Location: LCCOMB_X47_Y48_N12
\inst8|inst1|Mux1~5\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~5_combout\ = (\inst8|inst2|inst2~q\ & (!\inst8|inst1|Mux1~4_combout\)) # (!\inst8|inst2|inst2~q\ & ((!\inst8|inst2|inst~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000101001011111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst2~q\,
	datac => \inst8|inst1|Mux1~4_combout\,
	datad => \inst8|inst2|inst~q\,
	combout => \inst8|inst1|Mux1~5_combout\);

-- Location: LCCOMB_X47_Y48_N10
\inst8|inst1|Mux1~6\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~6_combout\ = (!\inst8|inst2|inst4~q\ & (!\inst8|inst2|inst3~q\ & \inst8|inst1|Mux1~5_combout\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000001100000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst2|inst4~q\,
	datac => \inst8|inst2|inst3~q\,
	datad => \inst8|inst1|Mux1~5_combout\,
	combout => \inst8|inst1|Mux1~6_combout\);

-- Location: LCCOMB_X46_Y48_N22
\inst8|inst1|Mux1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux1~combout\ = LCELL((\inst8|inst2|inst5~q\ & ((\inst8|inst1|Mux1~6_combout\))) # (!\inst8|inst2|inst5~q\ & (!\inst8|inst1|Mux1~3_combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000001010101",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst1|Mux1~3_combout\,
	datac => \inst8|inst1|Mux1~6_combout\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst8|inst1|Mux1~combout\);

-- Location: CLKCTRL_G5
\inst8|inst1|Mux1~clkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \inst8|inst1|Mux1~clkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \inst8|inst1|Mux1~clkctrl_outclk\);

-- Location: IOIBUF_X51_Y54_N22
\SW1~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_SW1,
	o => \SW1~input_o\);

-- Location: FF_X45_Y51_N27
\inst8|inst4|state.digit3\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \test_Sync~inputclkctrl_outclk\,
	asdata => \inst8|inst4|state.digit2~q\,
	sload => VCC,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst4|state.digit3~q\);

-- Location: LCCOMB_X45_Y51_N2
\inst8|inst4|state.digit0~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst4|state.digit0~0_combout\ = !\inst8|inst4|state.digit3~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111100001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst8|inst4|state.digit3~q\,
	combout => \inst8|inst4|state.digit0~0_combout\);

-- Location: FF_X45_Y51_N3
\inst8|inst4|state.digit0\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \test_Sync~inputclkctrl_outclk\,
	d => \inst8|inst4|state.digit0~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst4|state.digit0~q\);

-- Location: LCCOMB_X46_Y51_N12
\inst8|inst4|state.digit1~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst4|state.digit1~0_combout\ = !\inst8|inst4|state.digit0~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000011111111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datad => \inst8|inst4|state.digit0~q\,
	combout => \inst8|inst4|state.digit1~0_combout\);

-- Location: FF_X46_Y51_N13
\inst8|inst4|state.digit1\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \test_Sync~inputclkctrl_outclk\,
	d => \inst8|inst4|state.digit1~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst4|state.digit1~q\);

-- Location: FF_X45_Y51_N5
\inst8|inst4|state.digit2\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \test_Sync~inputclkctrl_outclk\,
	asdata => \inst8|inst4|state.digit1~q\,
	sload => VCC,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst8|inst4|state.digit2~q\);

-- Location: IOIBUF_X51_Y54_N29
\SW0~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_SW0,
	o => \SW0~input_o\);

-- Location: LCCOMB_X50_Y51_N0
\inst8|inst4|h~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst4|h~0_combout\ = (\SW1~input_o\) # ((\SW0~input_o\ & (!\inst8|inst4|state.digit2~q\)) # (!\SW0~input_o\ & ((!\inst8|inst4|state.digit1~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1011101110101111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \SW1~input_o\,
	datab => \inst8|inst4|state.digit2~q\,
	datac => \inst8|inst4|state.digit1~q\,
	datad => \SW0~input_o\,
	combout => \inst8|inst4|h~0_combout\);

-- Location: LCCOMB_X46_Y48_N2
\inst8|inst1|Mux2~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux2~0_combout\ = (!\inst8|inst2|inst3~q\ & ((\inst8|inst2|inst1~q\ & (!\inst8|inst2|inst2~q\)) # (!\inst8|inst2|inst1~q\ & (\inst8|inst2|inst2~q\ & !\inst8|inst2|inst~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000100110",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst1~q\,
	datab => \inst8|inst2|inst2~q\,
	datac => \inst8|inst2|inst~q\,
	datad => \inst8|inst2|inst3~q\,
	combout => \inst8|inst1|Mux2~0_combout\);

-- Location: LCCOMB_X46_Y48_N24
\inst8|inst1|Mux2~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst8|inst1|Mux2~1_combout\ = (!\inst8|inst2|inst4~q\ & (\inst8|inst1|Mux2~0_combout\ & !\inst8|inst2|inst5~q\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000001010000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst4~q\,
	datac => \inst8|inst1|Mux2~0_combout\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst8|inst1|Mux2~1_combout\);

-- Location: IOIBUF_X46_Y54_N15
\test_bcd_CLK~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_test_bcd_CLK,
	o => \test_bcd_CLK~input_o\);

-- Location: LCCOMB_X45_Y52_N22
\inst5|inst3|7~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst3|7~0_combout\ = !\inst5|inst3|7~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111100001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst3|7~q\,
	combout => \inst5|inst3|7~0_combout\);

-- Location: LCCOMB_X45_Y52_N24
\inst5|inst3|7~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst3|7~feeder_combout\ = \inst5|inst3|7~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst3|7~0_combout\,
	combout => \inst5|inst3|7~feeder_combout\);

-- Location: FF_X45_Y52_N25
\inst5|inst3|7\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \ALT_INV_test_bcd_CLK~input_o\,
	d => \inst5|inst3|7~feeder_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst3|7~q\);

-- Location: CLKCTRL_G13
\inst5|inst3|7~clkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \inst5|inst3|7~clkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \inst5|inst3|7~clkctrl_outclk\);

-- Location: LCCOMB_X45_Y51_N16
\inst5|inst3|31\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst3|31~combout\ = (!\inst5|inst3|19~q\ & !\inst5|inst3|11~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst3|19~q\,
	datad => \inst5|inst3|11~q\,
	combout => \inst5|inst3|31~combout\);

-- Location: FF_X45_Y51_N25
\inst5|inst3|11\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst3|ALT_INV_7~clkctrl_outclk\,
	asdata => \inst5|inst3|31~combout\,
	sload => VCC,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst3|11~q\);

-- Location: LCCOMB_X45_Y51_N18
\inst5|inst3|14~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst3|14~0_combout\ = \inst5|inst3|14~q\ $ (\inst5|inst3|11~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111111110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst3|14~q\,
	datad => \inst5|inst3|11~q\,
	combout => \inst5|inst3|14~0_combout\);

-- Location: FF_X45_Y51_N19
\inst5|inst3|14\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst3|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst3|14~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst3|14~q\);

-- Location: LCCOMB_X45_Y51_N0
\inst5|inst3|30\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst3|30~combout\ = (\inst5|inst3|14~q\ & \inst5|inst3|11~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst3|14~q\,
	datad => \inst5|inst3|11~q\,
	combout => \inst5|inst3|30~combout\);

-- Location: FF_X45_Y51_N1
\inst5|inst3|19\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst3|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst3|30~combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst3|19~q\);

-- Location: LCCOMB_X45_Y53_N18
\inst5|inst4|7~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst4|7~0_combout\ = !\inst5|inst4|7~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0011001100110011",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst4|7~q\,
	combout => \inst5|inst4|7~0_combout\);

-- Location: LCCOMB_X45_Y53_N28
\inst5|inst4|7~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst4|7~feeder_combout\ = \inst5|inst4|7~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst4|7~0_combout\,
	combout => \inst5|inst4|7~feeder_combout\);

-- Location: FF_X45_Y53_N29
\inst5|inst4|7\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst3|ALT_INV_19~q\,
	d => \inst5|inst4|7~feeder_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst4|7~q\);

-- Location: CLKCTRL_G10
\inst5|inst4|7~clkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \inst5|inst4|7~clkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \inst5|inst4|7~clkctrl_outclk\);

-- Location: LCCOMB_X46_Y51_N18
\inst5|inst4|31\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst4|31~combout\ = (!\inst5|inst4|19~q\ & !\inst5|inst4|11~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000110011",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst4|19~q\,
	datad => \inst5|inst4|11~q\,
	combout => \inst5|inst4|31~combout\);

-- Location: LCCOMB_X46_Y51_N6
\inst5|inst4|11~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst4|11~feeder_combout\ = \inst5|inst4|31~combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst4|31~combout\,
	combout => \inst5|inst4|11~feeder_combout\);

-- Location: FF_X46_Y51_N7
\inst5|inst4|11\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst4|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst4|11~feeder_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst4|11~q\);

-- Location: LCCOMB_X45_Y53_N24
\inst5|inst4|14~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst4|14~0_combout\ = \inst5|inst4|11~q\ $ (\inst5|inst4|14~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0011110000111100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst4|11~q\,
	datac => \inst5|inst4|14~q\,
	combout => \inst5|inst4|14~0_combout\);

-- Location: FF_X45_Y53_N25
\inst5|inst4|14\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst4|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst4|14~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst4|14~q\);

-- Location: LCCOMB_X45_Y53_N30
\inst5|inst4|30\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst4|30~combout\ = (\inst5|inst4|11~q\ & \inst5|inst4|14~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst4|11~q\,
	datad => \inst5|inst4|14~q\,
	combout => \inst5|inst4|30~combout\);

-- Location: FF_X45_Y53_N31
\inst5|inst4|19\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst4|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst4|30~combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst4|19~q\);

-- Location: LCCOMB_X44_Y53_N28
\inst5|inst5|7~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst5|7~0_combout\ = !\inst5|inst5|7~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111100001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst5|7~q\,
	combout => \inst5|inst5|7~0_combout\);

-- Location: LCCOMB_X44_Y53_N0
\inst5|inst5|7~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst5|7~feeder_combout\ = \inst5|inst5|7~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000011110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst5|7~0_combout\,
	combout => \inst5|inst5|7~feeder_combout\);

-- Location: FF_X44_Y53_N1
\inst5|inst5|7\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst4|ALT_INV_19~q\,
	d => \inst5|inst5|7~feeder_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst5|7~q\);

-- Location: CLKCTRL_G11
\inst5|inst5|7~clkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \inst5|inst5|7~clkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \inst5|inst5|7~clkctrl_outclk\);

-- Location: LCCOMB_X45_Y52_N26
\inst5|inst5|31\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst5|31~combout\ = (!\inst5|inst5|11~q\ & !\inst5|inst5|19~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst5|11~q\,
	datad => \inst5|inst5|19~q\,
	combout => \inst5|inst5|31~combout\);

-- Location: LCCOMB_X45_Y52_N30
\inst5|inst5|11~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst5|11~feeder_combout\ = \inst5|inst5|31~combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000011110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst5|31~combout\,
	combout => \inst5|inst5|11~feeder_combout\);

-- Location: FF_X45_Y52_N31
\inst5|inst5|11\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst5|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst5|11~feeder_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst5|11~q\);

-- Location: LCCOMB_X45_Y52_N10
\inst5|inst5|14~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst5|14~0_combout\ = \inst5|inst5|11~q\ $ (\inst5|inst5|14~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0101101001011010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst5|11~q\,
	datac => \inst5|inst5|14~q\,
	combout => \inst5|inst5|14~0_combout\);

-- Location: FF_X45_Y52_N11
\inst5|inst5|14\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst5|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst5|14~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst5|14~q\);

-- Location: LCCOMB_X45_Y52_N28
\inst5|inst5|30\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst5|30~combout\ = (\inst5|inst5|11~q\ & \inst5|inst5|14~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst5|11~q\,
	datad => \inst5|inst5|14~q\,
	combout => \inst5|inst5|30~combout\);

-- Location: FF_X45_Y52_N29
\inst5|inst5|19\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst5|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst5|30~combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst5|19~q\);

-- Location: LCCOMB_X44_Y52_N28
\inst5|inst6|7~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst6|7~0_combout\ = !\inst5|inst6|7~q\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111100001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst6|7~q\,
	combout => \inst5|inst6|7~0_combout\);

-- Location: LCCOMB_X44_Y52_N8
\inst5|inst6|7~feeder\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst6|7~feeder_combout\ = \inst5|inst6|7~0_combout\

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst5|inst6|7~0_combout\,
	combout => \inst5|inst6|7~feeder_combout\);

-- Location: FF_X44_Y52_N9
\inst5|inst6|7\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst5|ALT_INV_19~q\,
	d => \inst5|inst6|7~feeder_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst6|7~q\);

-- Location: CLKCTRL_G14
\inst5|inst6|7~clkctrl\ : fiftyfivenm_clkctrl
-- pragma translate_off
GENERIC MAP (
	clock_type => "global clock",
	ena_register_mode => "none")
-- pragma translate_on
PORT MAP (
	inclk => \inst5|inst6|7~clkctrl_INCLK_bus\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	outclk => \inst5|inst6|7~clkctrl_outclk\);

-- Location: LCCOMB_X49_Y51_N0
\inst5|inst6|31\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst6|31~combout\ = (!\inst5|inst6|11~q\ & !\inst5|inst6|19~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000001111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst6|11~q\,
	datad => \inst5|inst6|19~q\,
	combout => \inst5|inst6|31~combout\);

-- Location: FF_X49_Y51_N1
\inst5|inst6|11\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst6|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst6|31~combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst6|11~q\);

-- Location: LCCOMB_X49_Y51_N30
\inst5|inst6|14~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst6|14~0_combout\ = \inst5|inst6|14~q\ $ (\inst5|inst6|11~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000111111110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst6|14~q\,
	datad => \inst5|inst6|11~q\,
	combout => \inst5|inst6|14~0_combout\);

-- Location: FF_X49_Y51_N31
\inst5|inst6|14\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst6|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst6|14~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst6|14~q\);

-- Location: LCCOMB_X49_Y51_N28
\inst5|inst6|30\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst6|30~combout\ = (\inst5|inst6|14~q\ & \inst5|inst6|11~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datac => \inst5|inst6|14~q\,
	datad => \inst5|inst6|11~q\,
	combout => \inst5|inst6|30~combout\);

-- Location: FF_X49_Y51_N29
\inst5|inst6|19\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst5|inst6|ALT_INV_7~clkctrl_outclk\,
	d => \inst5|inst6|30~combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst5|inst6|19~q\);

-- Location: LCCOMB_X45_Y51_N20
\inst5|inst13~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst13~1_combout\ = (\inst8|inst4|state.digit3~q\ & ((\inst5|inst6|19~q\) # ((\inst8|inst4|state.digit1~q\ & \inst5|inst4|19~q\)))) # (!\inst8|inst4|state.digit3~q\ & (((\inst8|inst4|state.digit1~q\ & \inst5|inst4|19~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111100010001000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst4|state.digit3~q\,
	datab => \inst5|inst6|19~q\,
	datac => \inst8|inst4|state.digit1~q\,
	datad => \inst5|inst4|19~q\,
	combout => \inst5|inst13~1_combout\);

-- Location: LCCOMB_X45_Y51_N26
\inst5|inst13~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst13~0_combout\ = (\inst5|inst5|19~q\ & ((\inst8|inst4|state.digit2~q\) # ((\inst5|inst3|19~q\ & !\inst8|inst4|state.digit0~q\)))) # (!\inst5|inst5|19~q\ & (\inst5|inst3|19~q\ & ((!\inst8|inst4|state.digit0~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1010000011101100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst5|19~q\,
	datab => \inst5|inst3|19~q\,
	datac => \inst8|inst4|state.digit2~q\,
	datad => \inst8|inst4|state.digit0~q\,
	combout => \inst5|inst13~0_combout\);

-- Location: LCCOMB_X45_Y51_N24
\inst5|inst13\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst13~combout\ = (\inst5|inst13~1_combout\) # (\inst5|inst13~0_combout\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111111110101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst13~1_combout\,
	datad => \inst5|inst13~0_combout\,
	combout => \inst5|inst13~combout\);

-- Location: LCCOMB_X45_Y51_N30
\inst5|inst1~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst1~0_combout\ = (\inst5|inst5|7~q\ & ((\inst8|inst4|state.digit2~q\) # ((\inst5|inst3|7~q\ & !\inst8|inst4|state.digit0~q\)))) # (!\inst5|inst5|7~q\ & (((\inst5|inst3|7~q\ & !\inst8|inst4|state.digit0~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1000100011111000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst5|7~q\,
	datab => \inst8|inst4|state.digit2~q\,
	datac => \inst5|inst3|7~q\,
	datad => \inst8|inst4|state.digit0~q\,
	combout => \inst5|inst1~0_combout\);

-- Location: LCCOMB_X46_Y51_N26
\inst5|inst1~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst1~1_combout\ = (\inst8|inst4|state.digit3~q\ & ((\inst5|inst6|7~q\) # ((\inst5|inst4|7~q\ & \inst8|inst4|state.digit1~q\)))) # (!\inst8|inst4|state.digit3~q\ & (((\inst5|inst4|7~q\ & \inst8|inst4|state.digit1~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111100010001000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst4|state.digit3~q\,
	datab => \inst5|inst6|7~q\,
	datac => \inst5|inst4|7~q\,
	datad => \inst8|inst4|state.digit1~q\,
	combout => \inst5|inst1~1_combout\);

-- Location: LCCOMB_X46_Y51_N28
\inst5|inst1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst1~combout\ = (\inst5|inst1~0_combout\) # (\inst5|inst1~1_combout\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111111110101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst1~0_combout\,
	datad => \inst5|inst1~1_combout\,
	combout => \inst5|inst1~combout\);

-- Location: LCCOMB_X45_Y51_N4
\inst5|inst~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst~0_combout\ = (\inst5|inst3|14~q\ & (((\inst8|inst4|state.digit2~q\ & \inst5|inst5|14~q\)) # (!\inst8|inst4|state.digit0~q\))) # (!\inst5|inst3|14~q\ & (((\inst8|inst4|state.digit2~q\ & \inst5|inst5|14~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111001000100010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst3|14~q\,
	datab => \inst8|inst4|state.digit0~q\,
	datac => \inst8|inst4|state.digit2~q\,
	datad => \inst5|inst5|14~q\,
	combout => \inst5|inst~0_combout\);

-- Location: LCCOMB_X45_Y51_N6
\inst5|inst~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst~1_combout\ = (\inst5|inst4|14~q\ & ((\inst8|inst4|state.digit1~q\) # ((\inst5|inst6|14~q\ & \inst8|inst4|state.digit3~q\)))) # (!\inst5|inst4|14~q\ & (\inst5|inst6|14~q\ & (\inst8|inst4|state.digit3~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1110101011000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst4|14~q\,
	datab => \inst5|inst6|14~q\,
	datac => \inst8|inst4|state.digit3~q\,
	datad => \inst8|inst4|state.digit1~q\,
	combout => \inst5|inst~1_combout\);

-- Location: LCCOMB_X45_Y51_N22
\inst5|inst\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst~combout\ = (\inst5|inst~0_combout\) # (\inst5|inst~1_combout\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1110111011101110",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst~0_combout\,
	datab => \inst5|inst~1_combout\,
	combout => \inst5|inst~combout\);

-- Location: LCCOMB_X45_Y51_N8
\inst5|inst12~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst12~0_combout\ = (\inst8|inst4|state.digit0~q\ & (\inst8|inst4|state.digit2~q\ & ((\inst5|inst5|11~q\)))) # (!\inst8|inst4|state.digit0~q\ & ((\inst5|inst3|11~q\) # ((\inst8|inst4|state.digit2~q\ & \inst5|inst5|11~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1101110001010000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst4|state.digit0~q\,
	datab => \inst8|inst4|state.digit2~q\,
	datac => \inst5|inst3|11~q\,
	datad => \inst5|inst5|11~q\,
	combout => \inst5|inst12~0_combout\);

-- Location: LCCOMB_X46_Y51_N20
\inst5|inst12~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst12~1_combout\ = (\inst5|inst6|11~q\ & ((\inst8|inst4|state.digit3~q\) # ((\inst8|inst4|state.digit1~q\ & \inst5|inst4|11~q\)))) # (!\inst5|inst6|11~q\ & (\inst8|inst4|state.digit1~q\ & ((\inst5|inst4|11~q\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1110110010100000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst6|11~q\,
	datab => \inst8|inst4|state.digit1~q\,
	datac => \inst8|inst4|state.digit3~q\,
	datad => \inst5|inst4|11~q\,
	combout => \inst5|inst12~1_combout\);

-- Location: LCCOMB_X45_Y51_N10
\inst5|inst12\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst5|inst12~combout\ = (\inst5|inst12~0_combout\) # (\inst5|inst12~1_combout\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111111110101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst12~0_combout\,
	datad => \inst5|inst12~1_combout\,
	combout => \inst5|inst12~combout\);

-- Location: LCCOMB_X45_Y51_N12
\inst7|106~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|106~0_combout\ = (\inst5|inst1~combout\ & ((\inst5|inst12~combout\) # ((!\inst5|inst13~combout\ & !\inst5|inst~combout\)))) # (!\inst5|inst1~combout\ & (((!\inst5|inst~combout\ & \inst5|inst12~combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100111100000100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst13~combout\,
	datab => \inst5|inst1~combout\,
	datac => \inst5|inst~combout\,
	datad => \inst5|inst12~combout\,
	combout => \inst7|106~0_combout\);

-- Location: LCCOMB_X46_Y51_N30
\inst6|84~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst6|84~0_combout\ = (!\inst5|inst12~0_combout\ & (!\inst5|inst12~1_combout\ & ((\inst5|inst1~1_combout\) # (\inst5|inst1~0_combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000001010100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst12~0_combout\,
	datab => \inst5|inst1~1_combout\,
	datac => \inst5|inst1~0_combout\,
	datad => \inst5|inst12~1_combout\,
	combout => \inst6|84~0_combout\);

-- Location: LCCOMB_X46_Y51_N24
\inst6|84~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst6|84~1_combout\ = (\inst5|inst12~1_combout\ & (!\inst5|inst1~1_combout\ & (!\inst5|inst1~0_combout\))) # (!\inst5|inst12~1_combout\ & (\inst5|inst12~0_combout\ $ (((\inst5|inst1~1_combout\) # (\inst5|inst1~0_combout\)))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000001101010110",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst12~1_combout\,
	datab => \inst5|inst1~1_combout\,
	datac => \inst5|inst1~0_combout\,
	datad => \inst5|inst12~0_combout\,
	combout => \inst6|84~1_combout\);

-- Location: LCCOMB_X47_Y51_N22
\inst7|103~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|103~0_combout\ = (\inst5|inst~combout\ & ((!\inst6|84~1_combout\))) # (!\inst5|inst~combout\ & (\inst6|84~0_combout\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0011001110101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst6|84~0_combout\,
	datab => \inst6|84~1_combout\,
	datad => \inst5|inst~combout\,
	combout => \inst7|103~0_combout\);

-- Location: LCCOMB_X46_Y51_N14
\inst7|79~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|79~0_combout\ = (\inst5|inst~combout\ & (\inst5|inst13~combout\)) # (!\inst5|inst~combout\ & (((!\inst5|inst1~combout\ & \inst5|inst12~combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1000110110001000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst~combout\,
	datab => \inst5|inst13~combout\,
	datac => \inst5|inst1~combout\,
	datad => \inst5|inst12~combout\,
	combout => \inst7|79~0_combout\);

-- Location: LCCOMB_X46_Y51_N8
\inst7|78~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|78~0_combout\ = (\inst5|inst12~combout\ & ((\inst5|inst13~combout\) # ((\inst5|inst~combout\ & !\inst5|inst1~combout\)))) # (!\inst5|inst12~combout\ & (\inst5|inst~combout\ & (\inst5|inst1~combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111100000101000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst~combout\,
	datab => \inst5|inst1~combout\,
	datac => \inst5|inst12~combout\,
	datad => \inst5|inst13~combout\,
	combout => \inst7|78~0_combout\);

-- Location: LCCOMB_X47_Y51_N10
\inst7|17\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|17~combout\ = (((\inst8|inst2|inst4~q\) # (\inst8|inst2|inst5~q\)) # (!\inst8|inst1|Mux2~0_combout\)) # (!\inst8|inst4|state.digit0~q\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111111111110111",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst4|state.digit0~q\,
	datab => \inst8|inst1|Mux2~0_combout\,
	datac => \inst8|inst2|inst4~q\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst7|17~combout\);

-- Location: FF_X47_Y51_N11
\inst7|10\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|17~combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|10~q\);

-- Location: LCCOMB_X47_Y51_N2
\inst7|22~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|22~0_combout\ = (\inst8|inst1|Mux2~1_combout\ & (\inst8|inst4|state.digit1~q\)) # (!\inst8|inst1|Mux2~1_combout\ & ((\inst7|10~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110011110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst4|state.digit1~q\,
	datac => \inst7|10~q\,
	datad => \inst8|inst1|Mux2~1_combout\,
	combout => \inst7|22~0_combout\);

-- Location: FF_X47_Y51_N3
\inst7|23\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|22~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|23~q\);

-- Location: LCCOMB_X47_Y51_N16
\inst7|25~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|25~0_combout\ = (\inst8|inst1|Mux2~1_combout\ & (\inst8|inst4|state.digit2~q\)) # (!\inst8|inst1|Mux2~1_combout\ & ((\inst7|23~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100111111000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst4|state.digit2~q\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst7|23~q\,
	combout => \inst7|25~0_combout\);

-- Location: FF_X47_Y51_N17
\inst7|24\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|25~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|24~q\);

-- Location: LCCOMB_X47_Y51_N0
\inst7|45~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|45~0_combout\ = (\inst8|inst1|Mux2~1_combout\ & ((\inst8|inst4|state.digit3~q\))) # (!\inst8|inst1|Mux2~1_combout\ & (\inst7|24~q\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111110000001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst7|24~q\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst8|inst4|state.digit3~q\,
	combout => \inst7|45~0_combout\);

-- Location: FF_X47_Y51_N1
\inst7|47\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|45~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|47~q\);

-- Location: LCCOMB_X47_Y51_N26
\inst7|44~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|44~0_combout\ = (\inst7|47~q\ & (((\inst8|inst2|inst4~q\) # (\inst8|inst2|inst5~q\)) # (!\inst8|inst1|Mux2~0_combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1010101010100010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst7|47~q\,
	datab => \inst8|inst1|Mux2~0_combout\,
	datac => \inst8|inst2|inst4~q\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst7|44~0_combout\);

-- Location: FF_X47_Y51_N27
\inst7|46\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|44~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|46~q\);

-- Location: LCCOMB_X47_Y51_N8
\inst7|50~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|50~0_combout\ = (\inst7|46~q\ & ((\inst8|inst2|inst4~q\) # ((\inst8|inst2|inst5~q\) # (!\inst8|inst1|Mux2~0_combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000010110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst4~q\,
	datab => \inst8|inst1|Mux2~0_combout\,
	datac => \inst7|46~q\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst7|50~0_combout\);

-- Location: FF_X47_Y51_N9
\inst7|48\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|50~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|48~q\);

-- Location: LCCOMB_X47_Y51_N28
\inst7|51~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|51~0_combout\ = (\inst7|48~q\ & ((\inst8|inst2|inst4~q\) # ((\inst8|inst2|inst5~q\) # (!\inst8|inst1|Mux2~0_combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000010110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst4~q\,
	datab => \inst8|inst1|Mux2~0_combout\,
	datac => \inst7|48~q\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst7|51~0_combout\);

-- Location: FF_X47_Y51_N29
\inst7|49\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|51~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|49~q\);

-- Location: LCCOMB_X47_Y51_N14
\inst7|73~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|73~0_combout\ = (\inst7|49~q\ & ((\inst8|inst2|inst4~q\) # ((\inst8|inst2|inst5~q\) # (!\inst8|inst1|Mux2~0_combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000010110000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst4~q\,
	datab => \inst8|inst1|Mux2~0_combout\,
	datac => \inst7|49~q\,
	datad => \inst8|inst2|inst5~q\,
	combout => \inst7|73~0_combout\);

-- Location: FF_X47_Y51_N15
\inst7|75\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|73~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|75~q\);

-- Location: LCCOMB_X45_Y51_N28
\inst7|72~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|72~0_combout\ = (\inst5|inst~combout\ & (((\inst5|inst13~combout\ & \inst5|inst12~combout\)) # (!\inst5|inst1~combout\))) # (!\inst5|inst~combout\ & ((\inst5|inst13~combout\ & ((\inst5|inst12~combout\))) # (!\inst5|inst13~combout\ & 
-- (\inst5|inst1~combout\ & !\inst5|inst12~combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111001000100110",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst~combout\,
	datab => \inst5|inst1~combout\,
	datac => \inst5|inst13~combout\,
	datad => \inst5|inst12~combout\,
	combout => \inst7|72~0_combout\);

-- Location: LCCOMB_X47_Y51_N4
\inst7|72~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|72~1_combout\ = (\inst8|inst1|Mux2~1_combout\ & ((\inst7|72~0_combout\))) # (!\inst8|inst1|Mux2~1_combout\ & (\inst7|75~q\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111101000001010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst7|75~q\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst7|72~0_combout\,
	combout => \inst7|72~1_combout\);

-- Location: FF_X47_Y51_N5
\inst7|74\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|72~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|74~q\);

-- Location: LCCOMB_X47_Y51_N6
\inst7|78~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|78~1_combout\ = (\inst8|inst1|Mux2~1_combout\ & (\inst7|78~0_combout\)) # (!\inst8|inst1|Mux2~1_combout\ & ((\inst7|74~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100111111000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst7|78~0_combout\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst7|74~q\,
	combout => \inst7|78~1_combout\);

-- Location: FF_X47_Y51_N7
\inst7|76\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|78~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|76~q\);

-- Location: LCCOMB_X47_Y51_N30
\inst7|79~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|79~1_combout\ = (\inst8|inst1|Mux2~1_combout\ & (\inst7|79~0_combout\)) # (!\inst8|inst1|Mux2~1_combout\ & ((\inst7|76~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1010111110100000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst7|79~0_combout\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst7|76~q\,
	combout => \inst7|79~1_combout\);

-- Location: FF_X47_Y51_N31
\inst7|77\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|79~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|77~q\);

-- Location: FF_X47_Y51_N23
\inst7|103\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|103~0_combout\,
	asdata => \inst7|77~q\,
	sload => \inst8|inst1|ALT_INV_Mux2~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|103~q\);

-- Location: LCCOMB_X46_Y51_N0
\inst7|100~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|100~0_combout\ = (\inst5|inst1~combout\) # ((\inst5|inst~combout\ & (!\inst5|inst12~0_combout\ & !\inst5|inst12~1_combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111111100000010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst~combout\,
	datab => \inst5|inst12~0_combout\,
	datac => \inst5|inst12~1_combout\,
	datad => \inst5|inst1~combout\,
	combout => \inst7|100~0_combout\);

-- Location: LCCOMB_X47_Y51_N20
\inst7|100~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|100~1_combout\ = (\inst8|inst1|Mux2~1_combout\ & ((\inst7|100~0_combout\))) # (!\inst8|inst1|Mux2~1_combout\ & (\inst7|103~q\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000010101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst7|103~q\,
	datac => \inst7|100~0_combout\,
	datad => \inst8|inst1|Mux2~1_combout\,
	combout => \inst7|100~1_combout\);

-- Location: FF_X47_Y51_N21
\inst7|102\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|100~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|102~q\);

-- Location: LCCOMB_X47_Y51_N24
\inst7|106~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|106~1_combout\ = (\inst8|inst1|Mux2~1_combout\ & (\inst7|106~0_combout\)) # (!\inst8|inst1|Mux2~1_combout\ & ((\inst7|102~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100111111000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst7|106~0_combout\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst7|102~q\,
	combout => \inst7|106~1_combout\);

-- Location: FF_X47_Y51_N25
\inst7|104\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|106~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|104~q\);

-- Location: LCCOMB_X45_Y51_N14
\inst7|107~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|107~0_combout\ = (\inst5|inst~combout\ & (((\inst5|inst1~combout\ & \inst5|inst12~combout\)))) # (!\inst5|inst~combout\ & (!\inst5|inst13~combout\ & ((!\inst5|inst12~combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100000000010001",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst5|inst13~combout\,
	datab => \inst5|inst~combout\,
	datac => \inst5|inst1~combout\,
	datad => \inst5|inst12~combout\,
	combout => \inst7|107~0_combout\);

-- Location: LCCOMB_X47_Y51_N18
\inst7|107~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|107~1_combout\ = (\inst8|inst1|Mux2~1_combout\ & ((\inst7|107~0_combout\))) # (!\inst8|inst1|Mux2~1_combout\ & (\inst7|104~q\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111000010101010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst7|104~q\,
	datac => \inst7|107~0_combout\,
	datad => \inst8|inst1|Mux2~1_combout\,
	combout => \inst7|107~1_combout\);

-- Location: FF_X47_Y51_N19
\inst7|105\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|107~1_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|105~q\);

-- Location: LCCOMB_X47_Y51_N12
\inst7|123~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst7|123~0_combout\ = (\inst8|inst1|Mux2~1_combout\ & (\inst8|inst4|h~0_combout\)) # (!\inst8|inst1|Mux2~1_combout\ & ((\inst7|105~q\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100111111000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \inst8|inst4|h~0_combout\,
	datac => \inst8|inst1|Mux2~1_combout\,
	datad => \inst7|105~q\,
	combout => \inst7|123~0_combout\);

-- Location: FF_X47_Y51_N13
\inst7|124\ : dffeas
-- pragma translate_off
GENERIC MAP (
	is_wysiwyg => "true",
	power_up => "low")
-- pragma translate_on
PORT MAP (
	clk => \inst8|inst1|ALT_INV_Mux1~clkctrl_outclk\,
	d => \inst7|123~0_combout\,
	devclrn => ww_devclrn,
	devpor => ww_devpor,
	q => \inst7|124~q\);

-- Location: LCCOMB_X46_Y48_N20
inst101 : fiftyfivenm_lcell_comb
-- Equation(s):
-- \inst101~combout\ = (\inst8|inst1|Mux1~combout\ & ((\inst8|inst2|inst5~q\) # ((\inst8|inst2|inst4~q\) # (!\inst8|inst1|Mux2~0_combout\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100110010001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \inst8|inst2|inst5~q\,
	datab => \inst8|inst1|Mux1~combout\,
	datac => \inst8|inst1|Mux2~0_combout\,
	datad => \inst8|inst2|inst4~q\,
	combout => \inst101~combout\);

-- Location: IOIBUF_X49_Y54_N29
\KEY1~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_KEY1,
	o => \KEY1~input_o\);

-- Location: IOIBUF_X34_Y0_N29
\CLK_IN~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_CLK_IN,
	o => \CLK_IN~input_o\);

-- Location: IOIBUF_X46_Y54_N29
\KEY0~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_KEY0,
	o => \KEY0~input_o\);

-- Location: UNVM_X0_Y40_N40
\~QUARTUS_CREATED_UNVM~\ : fiftyfivenm_unvm
-- pragma translate_off
GENERIC MAP (
	addr_range1_end_addr => -1,
	addr_range1_offset => -1,
	addr_range2_offset => -1,
	is_compressed_image => "false",
	is_dual_boot => "false",
	is_eram_skip => "false",
	max_ufm_valid_addr => -1,
	max_valid_addr => -1,
	min_ufm_valid_addr => -1,
	min_valid_addr => -1,
	part_name => "quartus_created_unvm",
	reserve_block => "true")
-- pragma translate_on
PORT MAP (
	nosc_ena => \~QUARTUS_CREATED_GND~I_combout\,
	xe_ye => \~QUARTUS_CREATED_GND~I_combout\,
	se => \~QUARTUS_CREATED_GND~I_combout\,
	busy => \~QUARTUS_CREATED_UNVM~~busy\);

-- Location: ADCBLOCK_X43_Y52_N0
\~QUARTUS_CREATED_ADC1~\ : fiftyfivenm_adcblock
-- pragma translate_off
GENERIC MAP (
	analog_input_pin_mask => 0,
	clkdiv => 1,
	device_partname_fivechar_prefix => "none",
	is_this_first_or_second_adc => 1,
	prescalar => 0,
	pwd => 1,
	refsel => 0,
	reserve_block => "true",
	testbits => 66,
	tsclkdiv => 1,
	tsclksel => 0)
-- pragma translate_on
PORT MAP (
	soc => \~QUARTUS_CREATED_GND~I_combout\,
	usr_pwd => VCC,
	tsen => \~QUARTUS_CREATED_GND~I_combout\,
	chsel => \~QUARTUS_CREATED_ADC1~_CHSEL_bus\,
	eoc => \~QUARTUS_CREATED_ADC1~~eoc\);

-- Location: ADCBLOCK_X43_Y51_N0
\~QUARTUS_CREATED_ADC2~\ : fiftyfivenm_adcblock
-- pragma translate_off
GENERIC MAP (
	analog_input_pin_mask => 0,
	clkdiv => 1,
	device_partname_fivechar_prefix => "none",
	is_this_first_or_second_adc => 2,
	prescalar => 0,
	pwd => 1,
	refsel => 0,
	reserve_block => "true",
	testbits => 66,
	tsclkdiv => 1,
	tsclksel => 0)
-- pragma translate_on
PORT MAP (
	soc => \~QUARTUS_CREATED_GND~I_combout\,
	usr_pwd => VCC,
	tsen => \~QUARTUS_CREATED_GND~I_combout\,
	chsel => \~QUARTUS_CREATED_ADC2~_CHSEL_bus\,
	eoc => \~QUARTUS_CREATED_ADC2~~eoc\);

ww_STROBE <= \STROBE~output_o\;

ww_Data <= \Data~output_o\;

ww_D0 <= \D0~output_o\;

ww_D1 <= \D1~output_o\;

ww_D2 <= \D2~output_o\;

ww_D3 <= \D3~output_o\;

ww_CLK_OUT <= \CLK_OUT~output_o\;

ww_LATCH <= \LATCH~output_o\;

ww_A <= \A~output_o\;

ww_B <= \B~output_o\;

ww_C <= \C~output_o\;

ww_D <= \D~output_o\;

ww_h <= \h~output_o\;
END structure;


